#!/usr/bin/env bash
# Complete an AutoSolvate conda environment after `conda env create`.

set -euo pipefail

if [[ -z "${CONDA_PREFIX:-}" ]]; then
    echo "Error: activate the AutoSolvate conda environment before running this script." >&2
    exit 1
fi

python_executable="${CONDA_PREFIX}/bin/python"
if [[ ! -x "${python_executable}" ]]; then
    echo "Error: Python executable not found in active environment: ${python_executable}" >&2
    exit 1
fi

plugin_dir="${CONDA_PREFIX}/lib/openbabel/3.1.0"
if [[ ! -d "${plugin_dir}" ]]; then
    echo "Error: expected Open Babel plugin directory not found: ${plugin_dir}" >&2
    exit 1
fi

# Make the setting persistent for subsequent `conda activate` calls.
activate_dir="${CONDA_PREFIX}/etc/conda/activate.d"
mkdir -p "${activate_dir}"
printf '%s\n' \
    '#!/bin/sh' \
    '# Set by AutoSolvate post-create setup for conda-forge Open Babel.' \
    "export BABEL_LIBDIR=\"${plugin_dir}\"" \
    > "${activate_dir}/autosolvate-openbabel.sh"

export BABEL_LIBDIR="${plugin_dir}"

"${python_executable}" -m pip install --no-deps "molsimplify==1.8.0"

"${python_executable}" - <<'PY'
from openbabel import pybel
import molSimplify

print("Open Babel plugins and molSimplify imported successfully.")
PY

echo "Post-create setup complete. Reactivate the environment before running AutoSolvate"
echo "so BABEL_LIBDIR activation hook is loaded."

