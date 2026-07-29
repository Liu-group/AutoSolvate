"""Validate the external runtime expected by an AutoSolvate conda environment."""

from __future__ import annotations

import os
import sys
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path


def _check_environment() -> list[str]:
    errors: list[str] = []
    prefix_value = os.environ.get("CONDA_PREFIX")
    if not prefix_value:
        errors.append("CONDA_PREFIX is not set; activate the AutoSolvate conda environment.")
        prefix = None
    else:
        prefix = Path(prefix_value)

    babel_value = os.environ.get("BABEL_LIBDIR")
    if not babel_value:
        errors.append("BABEL_LIBDIR is not set; reactivate the conda environment.")
    elif not Path(babel_value).is_dir():
        errors.append(f"BABEL_LIBDIR does not exist: {babel_value}")
    elif prefix and prefix not in Path(babel_value).parents:
        errors.append(f"BABEL_LIBDIR is outside the active environment: {babel_value}")

    try:
        from openbabel import pybel  # noqa: F401
    except Exception as exc:
        errors.append(f"Open Babel Python/plugin import failed: {exc}")

    try:
        installed_version = version("molSimplify")
    except PackageNotFoundError:
        errors.append("molSimplify is not installed.")
    else:
        if installed_version != "1.8.0":
            errors.append(
                f"molSimplify 1.8.0 is required; found {installed_version}."
            )

    return errors


def main() -> int:
    errors = _check_environment()
    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        return 1

    print("AutoSolvate environment checks passed.")
    print(f"  CONDA_PREFIX={os.environ['CONDA_PREFIX']}")
    print(f"  BABEL_LIBDIR={os.environ['BABEL_LIBDIR']}")
    print("  molSimplify=1.8.0")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
