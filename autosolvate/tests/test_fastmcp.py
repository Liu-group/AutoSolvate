from __future__ import annotations

import asyncio
import json

import pytest

fastmcp = pytest.importorskip("fastmcp")
from fastmcp import Client

from autosolvate.fastmcp_server import mcp


def _payload(result):
    if hasattr(result, "data") and isinstance(result.data, dict):
        return result.data
    if hasattr(result, "content") and result.content:
        return json.loads(result.content[0].text)
    raise RuntimeError("Unexpected FastMCP call_tool result shape")


def test_fastmcp_water_acn_draft_and_cli_session():
    async def _run():
        async with Client(mcp) as client:
            response = await client.call_tool(
                "autosolvate_draft_water_acn_molar_mixture",
                {
                    "cube_size": 30,
                    "acetonitrile_molar_ratio": 0.3,
                    "water_molar_ratio": 0.7,
                    "charge_method": "bcc",
                },
            )
        payload = json.loads(response.content[0].text)
        normalized = payload["normalized_config"]

        assert normalized["system_type"] == "mixture"
        assert normalized["charge_method"] == "bcc"
        assert normalized["cube_size"] == 30
        names = [item["name"] for item in normalized["solvents"]]
        assert "acetonitrile" in names
        assert "water" in names

        async with Client(mcp) as client:
            started = await client.call_tool(
                "autosolvate_cli_session_start",
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
                "autosolvate_cli_session_send",
                {
                    "session_id": session_id,
                    "user_input": "exit",
                    "read_timeout_sec": 2.0,
                },
            )
            sent_payload = _payload(sent)
            assert isinstance(sent_payload["output"], str)

            closed = await client.call_tool(
                "autosolvate_cli_session_close",
                {
                    "session_id": session_id,
                },
            )
            closed_payload = _payload(closed)
            assert closed_payload["closed"] is True

    asyncio.run(_run())