from __future__ import annotations

from pathlib import Path


def resource_path(package: str, relative_path: str) -> str:
    """Return an on-disk path for a packaged resource.

    This replaces deprecated pkg_resources.resource_filename(). It relies on
    importlib.resources (stdlib) and works for both editable installs and
    normal installs.

    Parameters
    ----------
    package:
        Package name, e.g. "autosolvate".
    relative_path:
        Path relative to the package root, e.g. "data/ch3cn/ch3cn.pdb".

    Returns
    -------
    str
        Filesystem path to the requested resource.

    Notes
    -----
    For current AutoSolvate usage, resources are expected to be present as
    real files on disk (not inside a zip). In the future, if zip-safe
    packaging is needed, this function can be extended to copy resources to a
    temporary directory.
    """

    try:
        from importlib.resources import files  # py>=3.9

        return str(files(package).joinpath(relative_path))
    except Exception:
        # Python 3.8 fallback: importlib_resources backport.
        from importlib_resources import files  # type: ignore

        return str(files(package).joinpath(relative_path))


def autosolvate_resource(relative_path: str) -> str:
    return resource_path("autosolvate", relative_path)
