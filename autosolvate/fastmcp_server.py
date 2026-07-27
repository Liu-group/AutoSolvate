from __future__ import annotations

import argparse
import io
import logging
import os
import pty
import selectors
import signal
import sys
import threading
import time
import uuid
from contextlib import contextmanager, redirect_stderr, redirect_stdout
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

from autosolvate.clustergen import startclustergen
from autosolvate.generatetrajs import startmd


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


def _mcp_tool(*args, **kwargs):
    def decorator(func):
        if mcp is None:
            return func
        return mcp.tool(*args, **kwargs)(func)

    return decorator

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


def _get_session(session_id: str) -> _PtySession:
    with _SESSIONS_LOCK:
        if session_id not in _SESSIONS:
            raise ValueError(f"Unknown session_id: {session_id}")
        return _SESSIONS[session_id]


def _collect_generated_files(before: set[str], after: set[str], working_dir: str) -> list[str]:
    return sorted(os.path.abspath(os.path.join(working_dir, name)) for name in (after - before))


def _run_legacy_cli(entrypoint, argument_list: list[str], working_dir: str) -> Dict[str, Any]:
    stdout_buffer = io.StringIO()
    stderr_buffer = io.StringIO()
    os.makedirs(working_dir, exist_ok=True)
    with _pushd(working_dir):
        before = set(os.listdir("."))
        with redirect_stdout(stdout_buffer), redirect_stderr(stderr_buffer):
            entrypoint(argument_list)
        after = set(os.listdir("."))
    return {
        "working_dir": os.path.abspath(working_dir),
        "argument_list": argument_list,
        "stdout": stdout_buffer.getvalue(),
        "stderr": stderr_buffer.getvalue(),
        "generated_files": _collect_generated_files(before, after, working_dir),
    }


def _collect_reported_xyz_paths(stdout_text: str, working_dir: str) -> list[str]:
    paths = []
    for line in stdout_text.splitlines():
        candidate = line.strip()
        if not candidate.endswith(".xyz"):
            continue
        candidate_path = candidate if os.path.isabs(candidate) else os.path.join(working_dir, candidate)
        if os.path.exists(candidate_path):
            paths.append(os.path.abspath(candidate_path))
    return sorted(set(paths))


@_mcp_tool()
def autosolvate_boxgen_start(
    agent_mode: bool = True,
    working_dir: Optional[str] = None,
    initial_read_timeout_sec: float = 2.0,
) -> Dict[str, Any]:
    """
    AutoSolvate is a software package that generates the necessary input files (AMBER prmtop & inpcrd) for running a classical MD simulation of a solvated/condensed phase system. It streamlines the process of using multiple tools (such as antechamber/parmchk/packmol/tleap etc.) by providing a single command-line interface. 
    
    This "autosolvate_boxgen_start" tool starts an interactive wizard session for generating these input files. This session is managed by the AutoSolvate MCP itself. You only need to answer the questions asked by the wizard, and the tool will handle the rest. 

    The session is stateful and allows you to iteratively provide input and receive output until you complete the process or choose to exit. The tool returns a unique session_id that you can use to send further input to the wizard or to close the session when you're done.

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


@_mcp_tool()
def autosolvate_boxgen_send(
    session_id: str,
    user_input: str,
    read_timeout_sec: float = 1.5,
    max_bytes: int = 64000,
) -> Dict[str, Any]:
    """Send a single line to the wizard session to answer the question asked by the wizard, and read the response (or next question) from the wizard. You can call this tool iteratively until you complete the process or choose to exit."""
    session = _get_session(session_id)
    session.send(user_input)
    output = session.read(timeout_sec=read_timeout_sec, max_bytes=max_bytes)
    return {
        "session_id": session_id,
        "output": output.text,
        "timed_out": output.timed_out,
        "alive": session.is_alive(),
    }


@_mcp_tool()
def autosolvate_boxgen_close(session_id: str) -> Dict[str, Any]:
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


def _build_mdrun_argument_list(
    filename: str,
    temperature: float,
    pressure: float,
    stepsmmmin: int,
    stepsmmheat: int,
    stepsmmnve: int,
    stepsmmnpt: int,
    stepsqmmmmin: int,
    stepsqmmmheat: int,
    stepsqmmmnve: int,
    stepsqmmmnvt: int,
    charge: int,
    spinmultiplicity: int,
    functional: str,
    srun_use: bool,
    pmemduse: bool,
    dryrun: bool,
    freeze_solute: bool,
) -> list[str]:
    argument_list = [
        "-f", filename,
        "-t", str(temperature),
        "-p", str(pressure),
        "-i", str(stepsmmmin),
        "-m", str(stepsmmheat),
        "-b", str(stepsmmnve),
        "-n", str(stepsmmnpt),
        "-l", str(stepsqmmmmin),
        "-o", str(stepsqmmmheat),
        "-v", str(stepsqmmmnve),
        "-s", str(stepsqmmmnvt),
        "-q", str(charge),
        "-u", str(spinmultiplicity),
        "-k", functional,
    ]
    if srun_use:
        argument_list.append("-r")
    if pmemduse:
        argument_list.append("-x")
    if dryrun:
        argument_list.append("-d")
    if freeze_solute:
        argument_list.append("-z")
    return argument_list


@_mcp_tool()
def autosolvate_mdrun(
    filename: str,
    working_dir: str = ".",
    temperature: float = 300.0,
    pressure: float = 1.0,
    stepsmmmin: int = 2000,
    stepsmmheat: int = 10000,
    stepsmmnve: int = 0,
    stepsmmnpt: int = 300000,
    stepsqmmmmin: int = 250,
    stepsqmmmheat: int = 1000,
    stepsqmmmnve: int = 0,
    stepsqmmmnvt: int = 10000,
    charge: int = 0,
    spinmultiplicity: int = 1,
    functional: str = "b3lyp",
    srun_use: bool = False,
    pmemduse: bool = False,
    dryrun: bool = False,
    freeze_solute: bool = False,
) -> Dict[str, Any]:
    """Run the AutoSolvate MD driver with structured arguments."""
    argument_list = _build_mdrun_argument_list(
        filename=filename,
        temperature=temperature,
        pressure=pressure,
        stepsmmmin=stepsmmmin,
        stepsmmheat=stepsmmheat,
        stepsmmnve=stepsmmnve,
        stepsmmnpt=stepsmmnpt,
        stepsqmmmmin=stepsqmmmmin,
        stepsqmmmheat=stepsqmmmheat,
        stepsqmmmnve=stepsqmmmnve,
        stepsqmmmnvt=stepsqmmmnvt,
        charge=charge,
        spinmultiplicity=spinmultiplicity,
        functional=functional,
        srun_use=srun_use,
        pmemduse=pmemduse,
        dryrun=dryrun,
        freeze_solute=freeze_solute,
    )
    result = _run_legacy_cli(startmd, argument_list=argument_list, working_dir=working_dir)
    result.update({
        "filename": filename,
        "dryrun": dryrun,
    })
    return result


@_mcp_tool()
def autosolvate_clustergen(
    filename: str,
    trajname: str,
    working_dir: str = ".",
    startframe: int = 0,
    interval: int = 100,
    size: float = 4.0,
    srun_use: bool = False,
    spherical: bool = False,
) -> Dict[str, Any]:
    """Run the AutoSolvate cluster extraction driver with structured arguments."""
    argument_list = [
        "-f", filename,
        "-t", trajname,
        "-a", str(startframe),
        "-i", str(interval),
        "-s", str(size),
    ]
    if srun_use:
        argument_list.append("-r")
    if spherical:
        argument_list.append("-p")
    result = _run_legacy_cli(startclustergen, argument_list=argument_list, working_dir=working_dir)
    result["generated_files"] = sorted(
        set(result["generated_files"]) | set(_collect_reported_xyz_paths(result["stdout"], working_dir))
    )
    result.update({
        "filename": filename,
        "trajname": trajname,
    })
    return result


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
