![AutoSolvate](autosolvate/GUI/images/logo.png)

[//]: # (Badges)
[![Documentation Status](https://readthedocs.org/projects/autosolvate/badge/?version=latest)](https://autosolvate.readthedocs.io/en/latest/?badge=latest)
[![Anaconda-Server Badge](https://anaconda.org/liugroupemory/autosolvate/badges/version.svg)](https://anaconda.org/liugroupemory/autosolvate)
[![Anaconda-Server Badge](https://anaconda.org/liugroupemory/autosolvate/badges/license.svg)](https://anaconda.org/liugroupemory/autosolvate)
[![Anaconda-Server Badge](https://anaconda.org/liugroupemory/autosolvate/badges/installer/conda.svg)](https://conda.anaconda.org/liugroupemory)
[![Anaconda-Server Badge](https://anaconda.org/liugroupemory/autosolvate/badges/platforms.svg)](https://anaconda.org/liugroupemory/autosolvate)


Automated workflow for setting up explicit solvent quantum chemistry calculations for molecules.

### MCP JSON helper
AutoSolvate can emit a ready-to-use MCP server config (for GitHub Copilot, Claude, etc.).

```bash
autosolvate mcp_json --full > mcp.json
```

For a long-running HTTP MCP server:

```bash
autosolvate mcp_json --full --transport streamable-http --host 127.0.0.1 --port 8000 > mcp.json
autosolvate mcp_server --transport streamable-http --host 127.0.0.1 --port 8000
```

FastMCP is optional. Install it if you want the MCP server:

```bash
conda install -c conda-forge fastmcp
```

Or with pip:

```bash
pip install fastmcp
```

Options:
- `--server-name NAME`: customize the server name (default: `autosolvate-local`).
- `--indent N`: JSON indentation (default: 2).
- `--transport MODE`: choose `stdio` or `streamable-http`.

### Online Documentation
https://autosolvate.readthedocs.io/en/latest/

### Copyright

Copyright (c) 2022, Liu Group at Department of Chemistry, Emory University


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.3.
