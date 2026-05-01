from __future__ import annotations

import asyncio
import json

import pytest

fastmcp = pytest.importorskip("fastmcp")
from fastmcp import Client

from autosolvate.fastmcp_server import mcp
from autosolvate.tests import helper_functions as hp


def _payload(result):
    if hasattr(result, "data") and isinstance(result.data, dict):
        return result.data
    if hasattr(result, "content") and result.content:
        return json.loads(result.content[0].text)
    raise RuntimeError("Unexpected FastMCP call_tool result shape")


def test_fastmcp_boxgen_session_mdrun_and_clustergen(tmp_path):
    async def _run():
        async with Client(mcp) as client:
            started = await client.call_tool(
                "autosolvate_boxgen_start",
                {
                    "agent_mode": True,
                    "initial_read_timeout_sec": 2.5,
                },
            )
            start_payload = _payload(started)
            session_id = start_payload["session_id"]

            assert session_id
            assert "Welcome to the AutoSolvate Interactive Input Generator" in start_payload["initial_output"]

            sent = await client.call_tool(
                "autosolvate_boxgen_send",
                {
                    "session_id": session_id,
                    "user_input": "exit",
                    "read_timeout_sec": 2.0,
                },
            )
            sent_payload = _payload(sent)
            assert isinstance(sent_payload["output"], str)

            closed = await client.call_tool(
                "autosolvate_boxgen_close",
                {
                    "session_id": session_id,
                },
            )
            closed_payload = _payload(closed)
            assert closed_payload["closed"] is True

            mdrun = await client.call_tool(
                "autosolvate_mdrun",
                {
                    "filename": "water_solvated",
                    "working_dir": str(tmp_path / "mdrun"),
                    "dryrun": True,
                    "stepsmmheat": 10,
                    "stepsmmnpt": 10,
                    "stepsqmmmheat": 10,
                    "stepsqmmmnvt": 10,
                },
            )
            mdrun_payload = _payload(mdrun)
            generated_names = {path.split("/")[-1] for path in mdrun_payload["generated_files"]}
            assert "runMM.sh" in generated_names
            assert "runQMMMM.sh" in generated_names

            clustergen = await client.call_tool(
                "autosolvate_clustergen",
                {
                    "filename": hp.get_input_dir("water_solvated.prmtop"),
                    "trajname": hp.get_input_dir("water_solvated-qmmmnvt.netcdf"),
                    "working_dir": str(tmp_path / "clustergen"),
                    "startframe": 0,
                    "interval": 100,
                    "size": 4.0,
                },
            )
            clustergen_payload = _payload(clustergen)
            cluster_files = {path.split("/")[-1] for path in clustergen_payload["generated_files"]}
            assert "water_solvated-cutoutn-0.xyz" in cluster_files

    asyncio.run(_run())