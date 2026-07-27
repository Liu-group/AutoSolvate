from __future__ import annotations

import json
import os
from typing import Any, Dict, List, Optional, Sequence

from autosolvate.Common import SOLVENT_CLOSENESS, SOLVENT_DENSITY, SOLVENT_MW, WORKING_DIR
from autosolvate.utils.env_detection import available_qm_software, detect_amber_tools, detect_amberhome, cpu_count
from autosolvate.utils.interactive_input_utils import (
    InputAbort,
    apply_config_edit,
    ask_value,
    ask_yes_no,
    basename_no_ext,
    parse_cube_size,
    parse_float,
    parse_int,
)
from autosolvate.utils.solvent_calculation import (
    adjust_cube_size_for_density,
    calculate_solvent_number,
    calculate_solvent_numbers_from_molar_portions,
    calculate_solvent_number_from_density,
    calculate_solvent_numbers_from_volume_portions,
    calculate_solvent_numbers_from_weight_portions,
    cube_size_to_volume_m3,
    estimate_density_g_cm3,
    estimate_solvent_density,
    estimate_system_density,
    scale_solvent_numbers,
)
from autosolvate.utils.tools import check_multicomponent, check_transition_metal_complex, summarize_pdb_fragments, list_metal_atoms
from autosolvate.prompts.prompts import *
from autosolvate.utils.tools import determine_mw_from_xyz, determine_mw_from_pdb

MAX_RETRIES = 5
TARGET_DENSITY_G_CM3 = 1.0  # heuristic target for sanity checks


def _print_header():
    print("Welcome to the AutoSolvate Interactive Input Generator")
    print(f"Working directory: {os.getcwd()}")
    print("Control words: 'skip' to use default, 'exit' to abort\n")


def _detect_amber():
    amberhome = detect_amberhome()
    if amberhome:
        print(f"Default AMBERHOME: {amberhome}")
    else:
        print("Warning: AMBERHOME is not set. AmberTools is required for MCPB.")
    tools = detect_amber_tools(amberhome)
    for name, path in tools.items():
        if path:
            print(f"Detected {name}: {path}")
    if amberhome and all(v for v in tools.values()):
        print(f"Successfully detected AmberTools installation at {amberhome}")
    if amberhome and not amberhome.endswith("/"):
        amberhome += "/"
    return amberhome, tools


def _ask_cube_size(config: dict):
    cube = ask_value(
        "Enter the simulation box size in Angstrom.\n- ONE number for cubic (e.g., 50)\n- THREE numbers comma-separated for orthorhombic (e.g., 50,60,70)\n> ",
        parser=parse_cube_size,
        validator=lambda v: True,
    )
    config["cube_size"] = cube


def _ask_system_type() -> str:
    prompt = prompt_system_type_choice
    choice = ask_value(prompt, parser=parse_int, validator=lambda v: v in (1, 2))
    return "solute_solvent" if choice == 1 else "mixture"


def _auto_detect_kind(xyzfile: str) -> str:
    try:
        if check_multicomponent(xyzfile):
            return "complex"
        if check_transition_metal_complex(xyzfile):
            return "transition_metal_complex"
    except Exception:
        return "molecule"
    return "molecule"


def _exists_or_skip(value: str, allow_blank: bool = False) -> bool:
    if allow_blank and (value == "" or value is None):
        return True
    return os.path.exists(value)


def _ask_fragment_charges_with_summary(summary: List[tuple]) -> Dict[str, int]:
    """Prompt fragment charges with a printed summary of copies and atoms."""
    frag_charge: Dict[str, int] = {}
    for resname, copies, atoms in summary:
        total_atoms = copies * atoms
        print(f"Fragment {resname}: {copies} copies, {atoms} atoms each (total {total_atoms})")
    for resname, _, _ in summary:
        chg = ask_value(
            f"Charge for fragment {resname} (integer, default 0)\n> ",
            parser=parse_int,
            validator=lambda v: isinstance(v, int),
            default=0,
            allow_skip=True,
        )
        frag_charge[resname] = chg
    return frag_charge


def _ask_solutes(system_type: str) -> List[dict]:
    if system_type == "mixture":
        return []
    n = ask_value(
        prompt_n_solute_species,
        parser=parse_int,
        validator=lambda v: v >= 1,
    )
    solutes: List[dict] = []
    has_centered = False
    for i in range(1, n + 1):
        solute: Dict[str, Any] = {}
        path = ask_value(
            f"Solute #{i}: structure file path (.xyz/.pdb)\n> ",
            parser=lambda x: x.strip(),
            validator=lambda v: len(v) > 0 and _exists_or_skip(v),
        )
        solute["xyzfile"] = path
        default_name = basename_no_ext(path)
        name = ask_value(
            f"Solute #{i}: short name (default: {default_name})\n> ",
            parser=lambda x: x.strip(),
            validator=lambda v: True,
            default=default_name,
            allow_skip=True,
        )
        solute["name"] = name

        detected = _auto_detect_kind(path)
        metals = []
        if detected == "transition_metal_complex":
            metals = list_metal_atoms(path)
            print(f"Detected metal centers in solute #{i}: " + ", ".join(metals))
        type_backmap = {"molecule": "regular molecule", "transition_metal_complex": "transition metal complex", "complex": "multi-fragment/ion-pair"}
        kind_choice = ask_value(
            f"{prompt_classify_solute.rstrip()} (auto: {type_backmap[detected]})\n",
            parser      = parse_int,
            validator   = lambda v: v in (1, 2, 3),
            default     = {"molecule": 1, "transition_metal_complex": 2, "complex": 3}.get(detected, 1),
            allow_skip  = True,
        )
        kind_map = {1: "molecule", 2: "transition_metal_complex", 3: "complex"}
        solute_type = kind_map.get(kind_choice, detected)
        solute["solute_type"] = solute_type

        # Type-specific prompts
        if solute_type == "molecule":
            charge = ask_value(
                f"Solute #{i}: total charge (integer, default 0)\n> ",
                parser=parse_int,
                validator=lambda v: True,
                default=0,
                allow_skip=True,
            )
            solute["charge"] = charge
            spin = ask_value(
                f"Solute #{i}: spin multiplicity (integer, default 1)\n> ",
                parser=parse_int,
                validator=lambda v: v >= 1,
                default=1,
                allow_skip=True,
            )
            solute["spinmult"] = spin
        elif solute_type == "transition_metal_complex":
            charge = ask_value(
                f"Solute #{i}: total charge (integer, default 0)\n> ",
                parser=parse_int,
                validator=lambda v: True,
                default=0,
                allow_skip=True,
            )
            solute["charge"] = charge
            spin = ask_value(
                f"Solute #{i}: spin multiplicity (integer, default 1)\n> ",
                parser=parse_int,
                validator=lambda v: v >= 1,
                default=1,
                allow_skip=True,
            )
            solute["spinmult"] = spin
        elif solute_type == "complex":
            # complex: skip total charge/spin prompts; default neutral singlet
            solute["charge"] = 0
            solute["spinmult"] = 1
            if path.lower().endswith(".pdb"):
                summary = summarize_pdb_fragments(path)
                if summary:
                    frag_charge = _ask_fragment_charges_with_summary(summary)
                    solute["fragment_charge"] = frag_charge
                else:
                    print("Warning: could not parse fragments; defaulting fragment charges to 0.")
            else:
                print("Detected multi-fragment from non-PDB; defaulting fragment charges to 0 and spin to 1.")

        # shared prompts that should follow fragment charge gathering
        num = ask_value(
            f"Solute #{i}: number of copies (integer >=1, default 1)\n> ",
            parser=parse_int,
            validator=lambda v: v >= 1,
            default=1,
            allow_skip=True,
        )
        solute["number"] = num

        if solute_type == "molecule":
            has_params = ask_yes_no("Do you have frcmod/mol2/lib/off/prep for this solute? (yes/no, default no)\n> ", default=False)
            if has_params:
                have_frcmod = ask_yes_no("Do you have a frcmod file? (yes/no, default no)\n> ", default=False)
                if have_frcmod:
                    frcmod_path = ask_value(
                        "Path to frcmod file (or skip to ignore)\n> ",
                        parser=lambda x: x.strip(),
                        validator=lambda v: _exists_or_skip(v, allow_blank=True),
                        allow_skip=True,
                    )
                    if frcmod_path:
                        solute["frcmod"] = frcmod_path
                    main_param = ask_value(
                        "Provide ONE of mol2/off/lib/prep (or skip)\n> ",
                        parser=lambda x: x.strip(),
                        validator=lambda v: _exists_or_skip(v, allow_blank=True),
                        allow_skip=True,
                    )
                    if main_param:
                        ext = os.path.splitext(main_param)[1].lower().lstrip(".")
                        solute[ext] = main_param
                else:
                    # no params supplied
                    pass

        elif solute_type == "transition_metal_complex":
            solute["total_charge"] = charge
            solute["metal_charge"] = ask_value(
                "Metal center charge (integer, e.g., 2 for Cu(II) or 3 for Fe(III))\n> ",
                parser=parse_int,
                validator=lambda v: True,
            )
            chargefile = ask_value(
                prompt_ligand_chargefile,
                parser=lambda x: x.strip(),
                validator=lambda v: _exists_or_skip(v, allow_blank=True),
                default="",
                allow_skip=True,
            )
            solute["chargefile"] = chargefile or ""
            cutoff = ask_value(
                "Cutoff distance for coordinating ligands (Angstrom, default 2.8)\n> ",
                parser=parse_float,
                validator=lambda v: v > 0,
                default=2.8,
                allow_skip=True,
            )
            solute["cutoff"] = cutoff

        elif solute_type == "complex":
            # fragment handling already done; nothing else here
            pass

        # Centering
        if num == 1 and not has_centered:
            center = ask_yes_no("Center this solute in the box? (yes/no)\n> ", default=False)
            solute["centered"] = center
            has_centered = has_centered or center
        else:
            solute["centered"] = False

        solutes.append(solute)
    return solutes


def _ask_qm_settings(has_tmc: bool, config: dict, amberhome_default: Optional[str]):
    if not has_tmc:
        return
    print("At least one transition metal complex detected. MCPB requires QM settings.")
    detected = available_qm_software()
    avail_str = "\n".join([f"{k}: {v}" for k, v in detected.items() if v])
    if avail_str:
        print("Detected QM software:\n" + avail_str)
    short = ask_value(
        prompt_qm_software_choice,
        parser=lambda x: x.strip().lower(),
        validator=lambda v: v in {"orca", "gau", "g09", "g03", "gms", ""},
        default="orca",
        allow_skip=True,
    ) or "orca"
    qmexe = ask_value(
        prompt_qm_executable_path,
        parser=lambda x: x.strip(),
        validator=lambda v: True,
        default=detected.get(short) or "",
        allow_skip=True,
    )
    basis = ask_value(
        "Basis set (default lanl2dz)\n> ",
        parser=lambda x: x.strip(),
        validator=lambda v: True,
        default="lanl2dz",
        allow_skip=True,
    )
    method = ask_value(
        "DFT functional/method (default b3lyp)\n> ",
        parser=lambda x: x.strip(),
        validator=lambda v: True,
        default="b3lyp",
        allow_skip=True,
    )
    nprocs = ask_value(
        f"Number of processors (<= {cpu_count()}, default 1)\n> ",
        parser=parse_int,
        validator=lambda v: v >= 1,
        default=1,
        allow_skip=True,
    )
    maxcore = None
    if short == "orca":
        maxcore = ask_value(
            "ORCA maxcore in MB (default 4096)\n> ",
            parser=parse_int,
            validator=lambda v: v >= 1,
            default=4096,
            allow_skip=True,
        )
    opt = ask_yes_no("Run geometry optimization? (yes/no, default yes)\n> ", default=True)
    amberhome = amberhome_default or detect_amberhome() or ""
    amberhome = ask_value(
        prompt_amberhome_path,
        parser=lambda x: x.strip(),
        validator=lambda v: True,
        default=amberhome,
        allow_skip=True,
    )
    config.update({
        "software": short,
        "QMexe": qmexe,
        "basisset": basis,
        "method": method,
        "nprocs": nprocs,
        "opt": opt,
        # "amberhome": amberhome,   # TODO: automcpb treats antechamber, tleap, paths as $AMBERHOME/<executable>. Need to modify it into $AMBERHOME/bin/<executable> first.
    })
    if maxcore is not None:
        config["maxcore"] = maxcore
    if "cutoff" not in config:
        config["cutoff"] = 2.8
    config.setdefault("fakecharge", False)
    config.setdefault("mode", "A")


def _preset_solvent_defaults(name: str) -> Dict[str, Any]:
    defaults = {}
    if name in SOLVENT_DENSITY:
        defaults["density"] = SOLVENT_DENSITY[name]
    if name in SOLVENT_MW:
        defaults["molecular_weight"] = SOLVENT_MW[name]
    defaults["residue_name"] = name[:3].upper()
    return defaults


def _ask_solvents(cube_size: Any) -> List[dict]:
    n = ask_value(
        prompt_ask_solvent,
        parser=parse_int,
        validator=lambda v: v >= 1 and v <= 10,
    )
    method = 1
    if n > 1:
        method = ask_value(
            prompt_composition_method + "\n> ",
            parser=parse_int,
            validator=lambda v: v in (1, 2, 3, 4),
        )
    else:
        method = ask_value(
            prompt_single_component_method + "\n> ",
            parser=parse_int,
            validator=lambda v: v in (1, 2),
            default=1,
            allow_skip=True,
        )
    solvents: List[dict] = []
    for j in range(1, n + 1):
        comp: Dict[str, Any] = {}
        choice = ask_value(
            prompt_solvent_choice.format(i=j),
            parser=lambda x: x.strip().lower(),
            validator=lambda v: len(v) > 0,
        )
        if choice != "custom":
            comp["name"] = choice
            comp.update(_preset_solvent_defaults(choice))
        else:
            comp["xyzfile"] = ask_value(
                f"Component #{j} structure file (.xyz/.pdb)\n> ",
                parser=lambda x: x.strip(),
                validator=lambda v: len(v) > 0 and _exists_or_skip(v),
            )
            comp["name"] = os.path.splitext(os.path.basename(comp["xyzfile"]))[0]
            comp["residue_name"] = ask_value(
                f"Component #{j} residue name (2-3 chars)\n> ",
                parser=lambda x: x.strip(),
                validator=lambda v: len(v) > 0,
            )
            # optional params
            have_frcmod = ask_yes_no("Do you have a frcmod file? (yes/no, default no)\n If no, the forcefield parameter will be generated automatically from GAFF.\n> ", default=False)
            if have_frcmod:
                comp["frcmod"] = ask_value(
                    "frcmod path (or skip)\n> ",
                    parser=lambda x: x.strip(),
                    validator=lambda v: _exists_or_skip(v, allow_blank=True),
                    allow_skip=True,
                )
                main_param = ask_value(
                    "Provide ONE of mol2/off/lib/prep (or skip)\n> ",
                    parser=lambda x: x.strip(),
                    validator=lambda v: _exists_or_skip(v, allow_blank=True),
                    allow_skip=True,
                )
                if main_param:
                    ext = os.path.splitext(main_param)[1].lower().lstrip(".")
                    comp[ext] = main_param

        if method == 1:
            comp["number"] = ask_value(
                f"Number of molecules for component #{j} (integer >=0)\n> ",
                parser=parse_int,
                validator=lambda v: v >= 0,
            )
        elif method in [2, 3, 4] and n == 1:
            default_density = comp.get("density")
            if "density" not in comp and comp.get("name") in SOLVENT_DENSITY:
                default_density = SOLVENT_DENSITY[comp["name"]] / 1000.0
            if not default_density:
                density_prompt = f"Enter target density g/cm^3 for component #{j} (unit: g/cm^3, not kg/m^3; e.g., {TARGET_DENSITY_G_CM3})\n> "
            else:
                density_prompt = f"Enter target density g/cm^3 for component #{j} (unit: g/cm^3, not kg/m^3; default {default_density:.3f})\n> "
            density = ask_value(
                density_prompt,
                parser=parse_float,
                validator=lambda v: v > 0,
                default=default_density,
                allow_skip=True,
            )
            comp["density"] = density
            default_mw = comp.get("molecular_weight")
            if not default_mw and comp.get("name") in SOLVENT_MW:
                default_mw = SOLVENT_MW[comp["name"]]
            if not default_mw:
                if comp["xyzfile"].lower().endswith(".xyz"):
                    default_mw = determine_mw_from_xyz(comp["xyzfile"])
                elif comp["xyzfile"].lower().endswith(".pdb"):
                    default_mw = determine_mw_from_pdb(comp["xyzfile"])
            if not default_mw:
                mw_prompt = f"Enter molecular weight g/mol for component #{j}\n> "
            else:
                mw_prompt = f"Molecular weight g/mol for component #{j} (default {default_mw:.2f})\n> "
            mw = ask_value(
                mw_prompt,
                parser=parse_float,
                validator=lambda v: v > 0,
                default=default_mw,
                allow_skip=True,
            )
            comp["molecular_weight"] = mw
            # compute number from density
            volume_m3 = cube_size_to_volume_m3(cube_size)
            comp["number"] = calculate_solvent_number_from_density(mw, density, volume_m3)
            # set all other ratios to 1.0
            comp["weight_ratio"] = 1.0
            comp["volume_ratio"] = 1.0
            comp["molar_ratio"] = 1.0
            print(f"Computed molecule count: {comp['number']}")
        else:
            # density
            if "density" in comp:
                default_density = comp["density"]
            else:
                default_density = None
            if default_density is None:
                density_prompt = f"Enter the density g/cm^3 for component #{j} (unit: g/cm^3, not kg/m^3; e.g., {TARGET_DENSITY_G_CM3})\n> "
            else:
                density_prompt = f"Enter the density g/cm^3 for component #{j} (unit: g/cm^3, not kg/m^3; default {default_density:.3f})\n> "
            density = ask_value(
                density_prompt,
                parser=parse_float,
                validator=lambda v: v > 0,
                default=default_density,
                allow_skip=True,
            )
            comp["density"] = density
            # mw
            default_mw = comp.get("molecular_weight")
            if not default_mw and "xyzfile" in comp and os.path.exists(comp["xyzfile"]):
                try:
                    default_mw = determine_mw_from_xyz(comp["xyzfile"])
                except Exception:
                    default_mw = None
            mw = ask_value(
                f"Molecular weight g/mol for component #{j} (default {default_mw:.2f})\n> ",
                parser=parse_float,
                validator=lambda v: v > 0,
                default=default_mw,
                allow_skip=True,
            )
            comp["molecular_weight"] = mw
            ratio = ask_value(
                f"Enter ratio value for component #{j} (float >0)\n> ",
                parser=parse_float,
                validator=lambda v: v > 0,
            )
            if method == 2:
                comp["weight_ratio"] = ratio
            elif method == 3:
                comp["volume_ratio"] = ratio
            elif method == 4:
                comp["molar_ratio"] = ratio
        solvents.append(comp)

    # compute numbers if ratio-based (only when multiple components use ratios)
    if n > 1 and method in (2, 3, 4):
        volume_m3 = cube_size_to_volume_m3(cube_size)
        densities = [float(s["density"]) for s in solvents]
        mws = [float(s["molecular_weight"]) for s in solvents]
        if method == 2:
            ratios = [float(s["weight_ratio"]) for s in solvents]
            total = sum(ratios)
            ratios = [r / total for r in ratios]
            numbers = calculate_solvent_numbers_from_weight_portions(mws, densities, ratios, volume_m3)
        elif method == 3:
            ratios = [float(s["volume_ratio"]) for s in solvents]
            total = sum(ratios)
            ratios = [r / total for r in ratios]
            numbers = calculate_solvent_numbers_from_volume_portions(mws, densities, ratios, volume_m3)
        else:
            ratios = [float(s["molar_ratio"]) for s in solvents]
            total = sum(ratios)
            ratios = [r / total for r in ratios]
            numbers = calculate_solvent_numbers_from_molar_portions(mws, densities, ratios, volume_m3)
        for s, num in zip(solvents, numbers):
            s["number"] = num
    else:
        # number given; optionally infer missing density/mw from presets when single solvent
        for s in solvents:
            if "molecular_weight" not in s and s.get("name") in SOLVENT_MW:
                s["molecular_weight"] = SOLVENT_MW[s["name"]]
            if "density" not in s and s.get("name") in SOLVENT_DENSITY:
                s["density"] = SOLVENT_DENSITY[s["name"]] / 1000.0

    return solvents


def _density_sanity_checks(config: dict):
    solutes = config.get("solutes", [])
    solvents = config.get("solvents", [])
    cube = config.get("cube_size")
    if not solvents or cube is None:
        return
    current_density = estimate_system_density(solutes, solvents, cube)
    print(f"Estimated overall density: {current_density:.3f} g/cm^3")
    solvent_density = estimate_solvent_density(solvents, cube)
    print(f"Estimated solvent-only density: {solvent_density:.3f} g/cm^3")
    if solvent_density <= 0:
        print(
            "Note: solvent-only density could not be estimated (missing/invalid solvent number or molecular weight). "
            "Skipping density sanity check."
        )
        return
    deviation = abs(current_density - solvent_density) / solvent_density
    if deviation > 0.05:
        print(f"Warning: overall density considering solutes deviates significantly from solvent-only density by {deviation*100:.2f}%. Consider adjustments to solvent density to simulate dilute solutions.")
        if current_density > solvent_density:
            print("The system appears to be much denser than the pure solvent, this may cause issues in packmol initial structure generation and MD equilibration.")
        choice = ask_value(
            prompt_density_adjustment_choice + "\n> ",
            parser=parse_int,
            validator=lambda v: v in (1, 2, 3),
            default=3,
            allow_skip=True,
        )
        if choice == 1:
            new_cube = adjust_cube_size_for_density(cube, current_density, solvent_density)
            print(f"Preview new cube_size: {new_cube}")
            apply = ask_yes_no("Apply this cube_size? (yes/no)\n> ", default=True)
            if apply:
                config["cube_size"] = new_cube
        elif choice == 2:
            new_numbers = scale_solvent_numbers(solvents, current_density, solvent_density)
            print("Preview scaled solvent counts:")
            for s, n in zip(solvents, new_numbers):
                print(f"  {s.get('name', 'solvent')}: {s.get('number')} -> {n}")
            apply = ask_yes_no("Apply scaled counts? (yes/no)\n> ", default=True)
            if apply:
                for s, n in zip(solvents, new_numbers):
                    s["number"] = n


def _ask_packmol_and_output(config: dict):
    solvents = config.get("solvents", []) or []
    default_closeness = 2.0
    if len(solvents) == 1:
        name = solvents[0].get("name")
        if name in SOLVENT_CLOSENESS:
            default_closeness = SOLVENT_CLOSENESS[name]
            print(f"Single preset solvent '{name}' detected; suggested Packmol closeness: {default_closeness} Å (will be used as default).")
    closeness = ask_value(
        prompt_packmol_closeness + "\n> ",
        parser=parse_float,
        validator=lambda v: v > 0,
        default=default_closeness,
        allow_skip=True,
    )
    folder = ask_value(
        f"Output folder (default {os.getcwd()})\n> ",
        parser=lambda x: x.strip(),
        validator=lambda v: True,
        default=os.getcwd(),
        allow_skip=True,
    )
    config["closeness"] = closeness
    config["folder"] = folder


def _edit_mode(config: dict):
    print("Entering edit mode. Example: config[\"cube_size\"] = 60")
    for _ in range(MAX_RETRIES):
        cmd = input("Edit command (or type 'done' to finish): ")
        if cmd.strip().lower() in {"done", "exit"}:
            break
        if apply_config_edit(config, cmd):
            print("Applied.")
            continue
        else:
            print("Could not parse command; try again.")

def _review_and_edit(config: dict, agent_mode: bool = False):
    if agent_mode:
        print("Final configuration:")
        print(json.dumps(config, indent=2))
        confirm = ask_value(
            prompt_final_confirmation_agent,
            parser=lambda x: x.strip().lower(),
            validator=lambda v: v in {"proceed", "write", "exit"},
        )
        return confirm
    
    confirm = "no"
    while confirm.lower() != "yes":
        print("Current configuration:")
        print(json.dumps(config, indent=2))
        confirm = ask_value(
            "Type 'yes' to proceed with this configuration, 'edit' to modify, 'write_only' to save the JSON without execution, or 'exit' to abort.\n> ",
            parser=lambda x: x.strip().lower(),
            validator=lambda v: v in {"yes", "edit", "write_only", "exit"},
        )
        if confirm == "edit":
            _edit_mode(config)
        elif confirm == "exit":
            raise InputAbort("User aborted during review.")
        elif confirm == "write_only":
            break
    return confirm
    
def execute(config: dict):
    # Placeholder for execution logic
    print("Starting AutoSolvate execution (milestone output, details in logs)...", flush=True)
    from autosolvate.multicomponent import startmulticomponent_fromdata
    try:
        startmulticomponent_fromdata(config)
    except Exception:
        print("Execution failed. Check log files for details.", flush=True)
        raise
    else:
        print("Execution finished. Outputs and logs are in the working directory.", flush=True)

def run_wizard(agent_mode: bool = False) -> Dict[str, Any]:
    config: Dict[str, Any] = {}
    _print_header()
    amberhome, _ = _detect_amber()
    _ask_cube_size(config)
    system_type = _ask_system_type()
    config["system_type"] = system_type
    config["charge_method"] = "bcc"
    solutes = _ask_solutes(system_type)
    config["solutes"] = solutes
    has_tmc = any(s.get("solute_type") == "transition_metal_complex" for s in solutes)
    _ask_qm_settings(has_tmc, config, amberhome)
    solvents = _ask_solvents(config.get("cube_size"))
    config["solvents"] = solvents
    _density_sanity_checks(config)
    _ask_packmol_and_output(config)
    confirm = _review_and_edit(config, agent_mode=agent_mode)
    out_dir = config.get("folder", os.getcwd())
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    with open(os.path.join(out_dir, "wizard_input.json"), "w") as fh:
        json.dump(config, fh, indent=2)
    print(f"Wrote configuration to {os.path.join(out_dir, 'wizard_input.json')}")
    if confirm == "yes" and not agent_mode or confirm == "proceed" and agent_mode:
        execute(config)
    return config


def main():
    try:
        config = run_wizard(agent_mode=False)
    except InputAbort:
        print("Aborted by user.")
        return

if __name__ == "__main__":
    main()
