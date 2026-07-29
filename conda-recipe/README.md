# AutoSolvate conda recipe

Build the package from the repository root:

```bash
conda build conda-recipe/molsimplify -c conda-forge --no-anaconda-upload
conda build conda-recipe -c local -c conda-forge --no-anaconda-upload
```

Install the locally built package:

```bash
conda install -c local -c conda-forge autosolvate
```

The first recipe packages the exact molSimplify 1.8.0 compatibility version
required by AutoSolvate. Unlike the source-install helper, it does not run pip
after conda has finished or leave files that conda cannot track.

The AutoSolvate package installs activation hooks that set `BABEL_LIBDIR` from
the active environment. Validate an installation with:

```bash
autosolvate-doctor
```

The optional GUI and MCP integrations are not part of the base package. The MCP
server can be enabled with `conda install -c conda-forge "fastmcp>=3,<4"`.

The GitHub Actions workflow builds the package on every pull request. Tagged
releases are uploaded when the `ANACONDA_API_TOKEN` repository secret is set.
