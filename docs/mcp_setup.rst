MCP Setup
=================

AutoSolvate can run an MCP server for agent integrations (GitHub Copilot, Claude, etc.).
This section explains how to generate a ready-to-use MCP config and add it to your editor.

Prerequisites (optional)
------------------------

FastMCP is optional. Install it only if you want to run the MCP server:

.. code-block:: bash

   conda install -c conda-forge fastmcp

Or with pip:

.. code-block:: bash

   pip install fastmcp

Generate MCP JSON
-----------------

From the environment where AutoSolvate is installed:

.. code-block:: bash

   autosolvate mcp_json --full > mcp.json

For a long-running HTTP MCP server (recommended for frameworks that do not want
the server lifecycle tied to a stdio pipe), generate an HTTP config instead:

.. code-block:: bash

  autosolvate mcp_json --full --transport streamable-http --host 127.0.0.1 --port 8000 > mcp.json

This command:

- Detects the active Python executable and environment prefix.
- Populates `AMBERHOME`, `PATH`, and `LD_LIBRARY_PATH` when available.
- Emits warnings to `stderr` if key executables are missing.
- For HTTP configs, prints the recommended `autosolvate mcp_server ...` startup command to `stderr`.

Use the config in VS Code
-------------------------

Paste the generated server block into your VS Code (or other agent) MCP config (usually
`.vscode/mcp.json`). Example:

.. code-block:: json

   {
     "servers": {
       "autosolvate-local": {
         "type": "stdio",
         "command": "/path/to/python",
         "args": ["-m", "autosolvate", "mcp_server", "--transport", "stdio"],
         "env": {
           "PYTHONUNBUFFERED": "1",
           "AMBERHOME": "/path/to/env",
           "PATH": "/path/to/env/bin:${env:PATH}",
           "LD_LIBRARY_PATH": "/path/to/env/lib:${env:LD_LIBRARY_PATH}"
         }
       }
     },
     "inputs": []
   }

HTTP configuration example:

.. code-block:: json

   {
     "servers": {
       "autosolvate-local": {
         "type": "http",
         "url": "http://127.0.0.1:8000/mcp"
       }
     },
     "inputs": []
   }

Run the HTTP server separately
------------------------------

If you choose HTTP transport, start the server in a terminal before connecting
your MCP client:

.. code-block:: bash

   autosolvate mcp_server --transport streamable-http --host 127.0.0.1 --port 8000

`streamable-http` is usually the best choice for persistent integrations such as
LangChain adapters, because the AutoSolvate process is no longer tied to a stdio
pipe owned by the caller.

Notes
-----

- If `AMBERHOME` is not detected, set it manually to your AmberTools prefix.
- Ensure `packmol`, `antechamber`, and `tleap` are available in the same environment.
- You can rename the server via:

  .. code-block:: bash

     autosolvate mcp_json --full --server-name autosolvate-local
- You can customize the HTTP endpoint path if needed:

  .. code-block:: bash

     autosolvate mcp_json --full --transport streamable-http --host 127.0.0.1 --port 8000 --path /mcp
