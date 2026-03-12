from __future__ import annotations

import argparse
import copy
import json
import logging
import os
import pty
import selectors
import signal
import sys
import threading
import time
import uuid
from contextlib import contextmanager
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Optional

try:
    from fastmcp import FastMCP
except Exception as exc:  # pragma: no cover - import guarded for optional dependency
    FastMCP = None
    _FASTMCP_IMPORT_ERROR = exc
else:
    _FASTMCP_IMPORT_ERROR = None

from autosolvate.multicomponent import startmulticomponent_fromdata
from autosolvate.utils.env_detection import detect_amber_tools, detect_amberhome
from autosolvate.utils.inputparser import InputParser


if FastMCP is not None:
    mcp = FastMCP(
        name="AutoSolvate FastMCP",
        instructions=(
            "Tools for drafting and validating AutoSolvate multicomponent inputs, "
            "with an optional execution entrypoint."
        ),
    )
else:
    mcp = None

_LOGGER = logging.getLogger("autosolvate.mcp")
if not _LOGGER.handlers:
    logging.basicConfig(stream=sys.stderr, level=logging.INFO)


@dataclass
class _ReadResult:
    text: str
    timed_out: bool


class _PtySession:
    def __init__(self, argv: list[str], cwd: Optional[str] = None) -> None:
        self._argv = argv
        self._cwd = cwd
        self._pid: Optional[int] = None
        self._master_fd: Optional[int] = None
        self._buffer = ""

    @property
    def pid(self) -> int:
        if self._pid is None:
            raise RuntimeError("Session not started")
        return self._pid

    def _prepare_env(self) -> Dict[str, str]:
        env = os.environ.copy()
        repo_root = str(Path(__file__).resolve().parent.parent)
        py_path = env.get("PYTHONPATH", "")
        if py_path:
            env["PYTHONPATH"] = f"{repo_root}:{py_path}"
        else:
            env["PYTHONPATH"] = repo_root
        return env

    def start(self) -> None:
        if self._pid is not None:
            raise RuntimeError("Session already started")

        pid, master_fd = pty.fork()
        if pid == 0:
            env = self._prepare_env()
            if self._cwd is not None:
                os.chdir(self._cwd)
            os.execvpe(self._argv[0], self._argv, env)
            raise AssertionError("execvpe failed")

        self._pid = pid
        self._master_fd = master_fd
        os.set_blocking(master_fd, False)

    def is_alive(self) -> bool:
        if self._pid is None:
            return False
        try:
            pid, _ = os.waitpid(self._pid, os.WNOHANG)
        except ChildProcessError:
            return False
        return pid == 0

    def close(self, sig: int = signal.SIGTERM, timeout_sec: float = 2.0) -> None:
        if self._pid is None:
            return

        try:
            os.kill(self._pid, sig)
        except ProcessLookupError:
            pass

        deadline = time.time() + timeout_sec
        while time.time() < deadline:
            if not self.is_alive():
                break
            time.sleep(0.05)

        if self.is_alive():
            try:
                os.kill(self._pid, signal.SIGKILL)
            except ProcessLookupError:
                pass

        if self._master_fd is not None:
            try:
                os.close(self._master_fd)
            except OSError:
                pass

        self._pid = None
        self._master_fd = None

    def send(self, line: str) -> None:
        if self._master_fd is None:
            raise RuntimeError("Session not started")
        if not line.endswith("\n"):
            line += "\n"
        os.write(self._master_fd, line.encode("utf-8", errors="replace"))

    def read(self, timeout_sec: float = 2.0, max_bytes: int = 64000) -> _ReadResult:
        if self._master_fd is None:
            raise RuntimeError("Session not started")

        sel = selectors.DefaultSelector()
        sel.register(self._master_fd, selectors.EVENT_READ)

        start = time.time()
        read_bytes = 0
        timed_out = False
        try:
            while True:
                if (time.time() - start) >= timeout_sec:
                    timed_out = True
                    break
                if read_bytes >= max_bytes:
                    break
                events = sel.select(timeout=0.2)
                if not events:
                    continue
                try:
                    chunk = os.read(self._master_fd, 4096)
                except BlockingIOError:
                    continue
                except OSError:
                    break
                if not chunk:
                    break
                read_bytes += len(chunk)
                self._buffer += chunk.decode("utf-8", errors="replace")
        finally:
            try:
                sel.unregister(self._master_fd)
            except Exception:
                pass
            sel.close()

        text = self._buffer
        self._buffer = ""
        return _ReadResult(text=text, timed_out=timed_out)


_SESSIONS: Dict[str, _PtySession] = {}
_SESSIONS_LOCK = threading.Lock()


@contextmanager
def _pushd(path: str):
    original = os.getcwd()
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    try:
        yield
    finally:
        os.chdir(original)


def _validate_ratio(value: float, key: str) -> None:
    if value < 0 or value > 1:
        raise ValueError(f"{key} must be in [0, 1]. Got {value}")


def _normalize_with_parser(config: Dict[str, Any]) -> Dict[str, Any]:
    parser = InputParser()
    parser.read_dict(copy.deepcopy(config))
    parser.parse()
    return parser.data


def _get_session(session_id: str) -> _PtySession:
    with _SESSIONS_LOCK:
        if session_id not in _SESSIONS:
            raise ValueError(f"Unknown session_id: {session_id}")
        return _SESSIONS[session_id]


@mcp.tool()
def autosolvate_check_environment() -> Dict[str, Any]:
    """Check AMBERHOME detection and required AutoSolvate executables.

    Returns tool paths, missing tools list, and a boolean status.
    """
    amberhome = detect_amberhome()
    tools = detect_amber_tools(amberhome)
    missing = [name for name, path in tools.items() if not path]
    return {
        "amberhome": amberhome,
        "tools": tools,
        "all_tools_found": len(missing) == 0,
        "missing_tools": missing,
    }


@mcp.tool()
def autosolvate_cli_session_start(
    agent_mode: bool = True,
    working_dir: Optional[str] = None,
    initial_read_timeout_sec: float = 2.0,
) -> Dict[str, Any]:
    """Start a PTY-backed interactive AutoSolvate wizard session.

    Returns a session_id and the initial output from the wizard.
    """
    argv = [sys.executable, "-m", "autosolvate", "boxgen_interactive"]
    if agent_mode:
        argv.append("--agent-mode")

    session = _PtySession(argv=argv, cwd=working_dir)
    session.start()
    first = session.read(timeout_sec=initial_read_timeout_sec)
    session_id = str(uuid.uuid4())
    with _SESSIONS_LOCK:
        _SESSIONS[session_id] = session

    return {
        "session_id": session_id,
        "pid": session.pid,
        "initial_output": first.text,
        "timed_out": first.timed_out,
        "agent_mode": agent_mode,
        "working_dir": os.path.abspath(working_dir or os.getcwd()),
    }


@mcp.tool()
def autosolvate_cli_session_send(
    session_id: str,
    user_input: str,
    read_timeout_sec: float = 1.5,
    max_bytes: int = 64000,
) -> Dict[str, Any]:
    """Send a single line to the wizard session and return new output."""
    session = _get_session(session_id)
    session.send(user_input)
    output = session.read(timeout_sec=read_timeout_sec, max_bytes=max_bytes)
    return {
        "session_id": session_id,
        "output": output.text,
        "timed_out": output.timed_out,
        "alive": session.is_alive(),
    }


@mcp.tool()
def autosolvate_cli_session_read(
    session_id: str,
    timeout_sec: float = 1.0,
    max_bytes: int = 64000,
) -> Dict[str, Any]:
    """Read pending output from the wizard without sending input."""
    session = _get_session(session_id)
    output = session.read(timeout_sec=timeout_sec, max_bytes=max_bytes)
    return {
        "session_id": session_id,
        "output": output.text,
        "timed_out": output.timed_out,
        "alive": session.is_alive(),
    }


@mcp.tool()
def autosolvate_cli_session_close(session_id: str) -> Dict[str, Any]:
    """Terminate a wizard session and release its PTY resources."""
    with _SESSIONS_LOCK:
        if session_id not in _SESSIONS:
            return {
                "session_id": session_id,
                "closed": False,
                "message": "session not found",
            }
        session = _SESSIONS.pop(session_id)
    session.close()
    return {
        "session_id": session_id,
        "closed": True,
    }


@mcp.tool()
def autosolvate_draft_water_acn_molar_mixture(
    cube_size: float = 30.0,
    acetonitrile_molar_ratio: float = 0.3,
    water_molar_ratio: float = 0.7,
    charge_method: str = "bcc",
    include_empty_solute: bool = True,
) -> Dict[str, Any]:
    """Draft and normalize a water/acetonitrile molar-ratio mixture config.

    Uses AutoSolvate defaults for densities and molecular weights.
    Returns both the draft input and the normalized config.
    """
    _validate_ratio(acetonitrile_molar_ratio, "acetonitrile_molar_ratio")
    _validate_ratio(water_molar_ratio, "water_molar_ratio")
    ratio_sum = acetonitrile_molar_ratio + water_molar_ratio
    if abs(ratio_sum - 1.0) > 1e-6:
        raise ValueError(
            f"The sum of molar ratios must be 1.0. Got {ratio_sum:.8f}"
        )

    config: Dict[str, Any] = {
        "cube_size": cube_size,
        "charge_method": charge_method,
        "solvents": [
            {
                "name": "acetonitrile",
                "molar_ratio": acetonitrile_molar_ratio,
            },
            {
                "name": "water",
                "molar_ratio": water_molar_ratio,
            },
        ],
    }
    if include_empty_solute:
        config["solute"] = {}

    normalized = _normalize_with_parser(config)
    return {
        "draft_config": config,
        "normalized_config": normalized,
        "message": "Draft and normalization completed.",
    }


@mcp.tool()
def autosolvate_normalize_config(
    config: Dict[str, Any],
    write_full_json: bool = False,
    output_json: str = "autosolvate_input_full.json",
) -> Dict[str, Any]:
    """Normalize a config using InputParser and optionally write JSON to disk."""
    normalized = _normalize_with_parser(config)
    written_file = None
    if write_full_json:
        with open(output_json, "w", encoding="utf-8") as handle:
            json.dump(normalized, handle, indent=2)
        written_file = os.path.abspath(output_json)
    return {
        "normalized_config": normalized,
        "written_file": written_file,
    }


@mcp.tool()
def autosolvate_run_multicomponent(
    config: Dict[str, Any],
    output_dir: str = ".",
    execute: bool = False,
) -> Dict[str, Any]:
    """Write draft + normalized JSON, optionally execute multicomponent build.

    If execute=False, only JSON files are written; no external tools are run.
    """
    normalized = _normalize_with_parser(config)
    with _pushd(output_dir):
        with open("wizard_input.json", "w", encoding="utf-8") as handle:
            json.dump(config, handle, indent=2)
        with open("autosolvate_input_full.json", "w", encoding="utf-8") as handle:
            json.dump(normalized, handle, indent=2)
        generated_files = [
            os.path.abspath("wizard_input.json"),
            os.path.abspath("autosolvate_input_full.json"),
        ]

        if execute:
            startmulticomponent_fromdata(config)
            for name in os.listdir("."):
                if name.endswith((".prmtop", ".inpcrd", ".pdb", ".mol2", ".frcmod")):
                    generated_files.append(os.path.abspath(name))

    return {
        "executed": execute,
        "output_dir": os.path.abspath(output_dir),
        "generated_files": sorted(set(generated_files)),
        "normalized_config": normalized,
    }


def main(argv: list[str] | None = None) -> None:
    if FastMCP is None:
        _LOGGER.error(
            "FastMCP is not installed. Install via conda-forge (recommended): 'conda install -c conda-forge fastmcp' "
            "or pip: 'pip install fastmcp'."
        )
        if _FASTMCP_IMPORT_ERROR is not None:
            _LOGGER.error("FastMCP import error: %s", _FASTMCP_IMPORT_ERROR)
        raise SystemExit(1)
    parser = argparse.ArgumentParser(description="AutoSolvate FastMCP server")
    parser.add_argument(
        "--transport",
        default="stdio",
        choices=["stdio", "streamable-http", "sse", "http"],
        help="FastMCP transport mode",
    )
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--port", type=int, default=8000)
    args = parser.parse_args(argv)

    _LOGGER.info("Starting AutoSolvate FastMCP server")
    _LOGGER.info("Transport=%s Host=%s Port=%s", args.transport, args.host, args.port)

    if args.transport == "stdio":
        try:
            mcp.run(transport="stdio", show_banner=False)
        except Exception:
            _LOGGER.exception("FastMCP stdio server failed")
            raise
        _LOGGER.error("FastMCP stdio server exited unexpectedly")
        return

    try:
        mcp.run(transport=args.transport, host=args.host, port=args.port, show_banner=False)
    except Exception:
        _LOGGER.exception("FastMCP server failed")
        raise
    _LOGGER.error("FastMCP server exited unexpectedly")


if __name__ == "__main__":
    main()
