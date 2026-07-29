# AutoSolvate conda recipe

Build the package from the repository root:

```bash
conda build conda-recipe -c conda-forge --no-anaconda-upload
```

Install the locally built package:

```bash
conda install -c local -c conda-forge autosolvate
```

The optional GUI and MCP integrations are not part of the base package. The MCP
server can be enabled with `conda install -c conda-forge "fastmcp>=3,<4"`.

The GitHub Actions workflow builds the package on every pull request. Tagged
releases are uploaded when the `ANACONDA_API_TOKEN` repository secret is set.
