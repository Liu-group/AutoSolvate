## feat(resp): Add NumPy vdwsurface Fallback, Unify `gen_esp`, and Relax QC Path Checks

- **Add** `vdwsurface_numpy.py` as an internal pure-Python fallback.
- **Update** `gen_esp.py` to:
    - Prefer `pyvdwsurface` when available.
    - Fallback to internal NumPy vdwsurface when unavailable.
    - Set default `ESP_GRID_DENSITY` to 3.
- **Remove** redundant `autosolvate/gen_esp.py` and consolidate usage to `resp_classes.gen_esp`.
- **Update** `FFmetalcomplex.py` import to `from .resp_classes import gen_esp`.
- **Relax** early executable validation in `resp_orca.py` and `resp_gamess.py`:
    - Warn during init instead of hard-fail.
    - Enforce executable checks only at actual execution time.
- **In ORCA flow**, skip `orca_vpot` when existing ESP output is already present.

_This improves test robustness in prebuilt-result scenarios and reduces RESP variability from overly dense ESP grids._