# Copilot instructions for AutoSolvate

AutoSolvate is a Python package/CLI/GUI for building explicitly solvated systems, running Amber MD/QM-MM workflows, and extracting microsolvated clusters.

## Big picture (where to look)
- **Entrypoint/CLI routing**: `autosolvate/__main__.py` (no args launches GUI; subcommands: `boxgen`, `boxgen_metal`, `boxgen_multicomponent`, `mdrun`, `clustergen`, `boxgen_interactive`).
- **Core “build solvated box” pipeline**: `autosolvate/autosolvate.py` (`solventBoxBuilder`, `AmberParamsBuilder`, and `startboxgen`).
- **MD + QM/MM driver**: `autosolvate/generatetrajs.py` (writes Amber input files; runs `sander`/`pmemd.cuda`; supports `--dryrun` and optional Slurm `srun`).
- **Cluster extraction**: `autosolvate/clustergen.py` (preferred path is `clustergen_pytraj` using `pytraj`/`parmed`).
- **Domain model objects**: `autosolvate/molecule/` (e.g., `SolventBox` and `AMBER_SOLVENTBOX_DICT` in `autosolvate/molecule/solventbox.py`).
- **External tool wrappers (“dockers”)**: `autosolvate/dockers/` (e.g., `Antechamber*Docker`, `PackmolDocker`, `TleapDocker`). These are thin wrappers around CLI programs.

## Key conventions and patterns (project-specific)
- **“Docker” classes are not containers**: they are wrappers around external executables and should:
  - implement `check_system()`, `predict_output()` (populate `self.output_files`), `generate_cmd()`, `execute()`, `check_output()`, `process_output()`, `run()`.
  - prefer `GeneralDocker.resolve_executable(..., amberhome)` when calling AmberTools (`antechamber`, `tleap`, `parmchk2`).
  - log to `autosolvate.log` (see `autosolvate/dockers/general_docker.py`).
- **Packmol path quirk**: `PackmolDocker` copies input PDBs into the workfolder because long paths can break packmol.
- **HPC/Slurm awareness**: many commands support a `srun` prefix (CLI options `-r/--srunuse`); do not assume `srun` exists (see docs note in `docs/installation.rst`).
- **Resource files** (solvent templates etc.) are accessed via `autosolvate.utils.resources.autosolvate_resource()`.
- **Interactive wizard**: `autosolvate/input_wizard.py` + `autosolvate/cli/boxgen_interactive.py` generates a JSON config; it detects AmberTools/QM executables (`autosolvate/utils/env_detection.py`).

## Developer workflows (verified in repo)
- **Recommended env**: `conda env create -f devtools/conda-envs/test_env.yaml` then `conda activate autosolvate` (see `docs/installation.rst`).
- **Install from source** (in the env): `python setup.py install`.
- **Run tests**: `pytest` (tests live in `autosolvate/tests/`; fixtures copy `autosolvate/tests/inputs` into tmpdir).
- **Sphinx docs**: sources in `docs/` (see `docs/index.rst`).

## External dependencies and integration points
- **AmberTools**: `antechamber`, `parmchk2`, `tleap`, `sander` are assumed available (often via `AMBERHOME`).
- **Packmol**: used for building custom solvent boxes.
- **OpenBabel**: used heavily for IO and fragment detection (`openbabel.pybel`).
- **Optional QM backends**: Gaussian/ORCA/GAMESS/TeraChem for RESP/QM-MM pathways (see `autosolvate/resp_classes/` and `autosolvate/dockers/terachem_docker.py`).

## Practical guidance for changes
- Prefer extending existing pipelines (CLI `start...` functions and docker wrappers) instead of introducing new execution layers.
- When adding a new external-tool wrapper, follow the existing docker interface and ensure `predict_output()` lists generated files for milestone printing.
- When changing CLI behavior, update `autosolvate/__main__.py` and keep the “no args launches GUI” behavior.
- Keep changes minimal: this codebase mixes older script-style modules (e.g., `generatetrajs.py`) with newer structured utilities (wizard + env detection).
