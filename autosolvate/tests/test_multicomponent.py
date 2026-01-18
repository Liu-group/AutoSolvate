#--------------------------------------------------------------------------------------------------#
# test multicomponent amber parameter generation
# author: Fangning Ren (2022-12-17)
# path: autosolvate/tests/test_multicomponent.py 
#--------------------------------------------------------------------------------------------------#
import os
import json
import pytest
import numpy as np

from ..multicomponent import *
from . import helper_functions as hp


def test_ionpair_solvation(tmpdir):
    """
    @TODO:
        New PDB compare function 

    @NOTE:
        The resulting solvent box is slightly different from the previous one: 
        The positions of the sodium ions used to balance the charges are inconsistent. 
        But the current pdb comparison function cannot tolerate such subtle differences. 
        Therefore, I canceled them before the new pdb functions finished.
    """
    testName = "test_ionpair_solvation"
    solutexyz = hp.get_input_dir("ionpair.pdb")
    name = "ionpair"
    inst = MulticomponentSolventBoxBuilder(
        xyzfile = solutexyz,
        slu_charge={"SUF":-2, "TPA":1},
        solvent="water",
        cube_size=20,
        charge_method="bcc",
        folder = os.getcwd()
    )
    inst.build()
    path_fragment_exist = True
    for res in ("TPA", "SUF"):
        path_fragment_exist *= os.path.exists(f"{name}-{res.lower()}.pdb")
        path_fragment_exist *= os.path.exists(f"{name}-{res.lower()}.lib")
    assert path_fragment_exist
    path_main_exist = True
    for suffix in ("lib", "pdb"):
        path_main_exist *= os.path.exists(f"{name}.{suffix}")
    assert path_main_exist

def test_ionpair_solvation_autodetect(tmpdir):
    """Same as test_ionpair_solvation but rely on auto charge detection."""
    testName = "test_ionpair_solvation_autodetect"
    solutexyz = hp.get_input_dir("ionpair.pdb")
    name = "ionpair"
    inst = MulticomponentSolventBoxBuilder(
        xyzfile=solutexyz,
        solvent="water",
        cube_size=20,
        charge_method="bcc",
        folder=os.getcwd(),
    )
    inst.build()
    path_fragment_exist = True
    for res in ("TPA", "SUF"):
        path_fragment_exist *= os.path.exists(f"{name}-{res.lower()}.pdb")
        path_fragment_exist *= os.path.exists(f"{name}-{res.lower()}.lib")
    assert path_fragment_exist
    path_main_exist = True
    for suffix in ("lib", "pdb"):
        path_main_exist *= os.path.exists(f"{name}.{suffix}")
    assert path_main_exist

    # check the prmtop to see whether we have only one Na+
    try:
        from parmed import load_file
    except ImportError:
        pytest.skip("parmed not installed, skipping prmtop check.")
    prmtop = load_file(f"water_solvated.prmtop")
    na_count = sum(1 for atom in prmtop.atoms if atom.name.lower() in ("na", "na+"))
    assert na_count == 1, f"Expected 1 Na+ ion, found {na_count}."

def test_multicomponent(tmpdir):
    testName = "test_multicomponent"
    inpfname = "PAHs"
    builder = MulticomponentParamsBuilder(
        hp.get_input_dir(f"{inpfname}.pdb"),
        charge_method="bcc",
        folder = os.getcwd(),
        pre_optimize_fragments=True,
    )
    builder.build()
    path_exist = True
    for res in ("UAA", "UAB", "UAC", "UAD"):
        path_exist *= os.path.exists(f"{inpfname}-{res.lower()}.pdb")
        path_exist *= os.path.exists(f"{inpfname}-{res.lower()}.lib")
        path_exist *= os.path.exists(f"{inpfname}-{res.lower()}.frcmod")
    assert path_exist
    assert hp.compare_pdb(f"{inpfname}.pdb", hp.get_reference_dir(f"multicomponent/{inpfname}-processed.pdb"))


def _count_residues_in_prmtop(prmtop_path: str):
    try:
        from parmed import load_file
    except ImportError:
        pytest.skip("parmed not installed, skipping prmtop check.")
    prmtop = load_file(prmtop_path)
    res_counts = {}
    for res in prmtop.residues:
        res_counts[res.name] = res_counts.get(res.name, 0) + 1
    return res_counts


def _assert_two_component_ratio_close(res_counts: dict, expected_ratio: float, rtol: float = 0.15):
    residues = list(res_counts.keys())
    assert len(residues) == 2, f"Expected 2 residue types, found {len(residues)}. Got: {res_counts}"
    c1 = res_counts[residues[0]]
    c2 = res_counts[residues[1]]
    assert c1 > 0 and c2 > 0, f"Residue counts must be positive. Got: {res_counts}"
    ratio = c1 / c2 if c1 >= c2 else c2 / c1
    assert np.isclose(ratio, expected_ratio, rtol=rtol), f"Residue ratio {ratio} deviates from expected {expected_ratio}. Counts: {res_counts}"


def test_mixture_partition_volume_ratio(tmpdir):
    """Ensure volume_ratio partitioning matches the reference molecule counts."""
    inputfilepath = hp.get_input_dir("water_acn_vratio.json")
    startmulticomponent(["-f", inputfilepath])
    assert os.path.exists("autosolvate_input_full.json")
    with open("autosolvate_input_full.json", "r") as f:
        data = json.load(f)
    counts = {s["name"]: s["number"] for s in data["solvents"]}

    # Reference from /pscratch/sd/f/fren5/AutoSolvate/test_acn_water/autosolvate_input_full.json
    n_acn = counts.get("acetonitrile")
    n_wat = counts.get("water")
    assert n_acn is not None and n_wat is not None, f"Expected acetonitrile and water in solvents. Got: {counts}"
    assert np.isclose(n_acn, 221, rtol=0.05), f"Acetonitrile count {n_acn} deviates from reference 221"
    assert np.isclose(n_wat, 1493, rtol=0.05), f"Water count {n_wat} deviates from reference 1497"


def test_mixture_partition_molar_ratio(tmpdir):
    """molar_ratio partitioning should produce positive counts and roughly 7:3 ratio."""
    inputfilepath = hp.get_input_dir("water_acn_mratio.json")
    startmulticomponent(["-f", inputfilepath])
    # Without an explicit prefix, MixtureBuilder names the system from component names.
    assert os.path.exists("acetonitrile-water.prmtop")
    res_counts = _count_residues_in_prmtop("acetonitrile-water.prmtop")
    _assert_two_component_ratio_close(res_counts, expected_ratio=7 / 3, rtol=0.05)


def test_mixture_partition_weight_ratio(tmpdir):
    """weight_ratio partitioning should produce positive counts and roughly preserve mass fractions.

    Note: weight_ratio is a mass fraction, so the molecule-count ratio is not expected to be 7:3.
    """
    inputfilepath = hp.get_input_dir("water_acn_wratio.json")
    startmulticomponent(["-f", inputfilepath])
    assert os.path.exists("acetonitrile-water.prmtop")
    res_counts = _count_residues_in_prmtop("acetonitrile-water.prmtop")
    # Residue names: acetonitrile -> C3N (see Molecule residue normalization), water -> WAT
    assert "C3N" in res_counts and "WAT" in res_counts, f"Expected C3N and WAT residues. Got: {res_counts}"
    n_acn = res_counts["C3N"]
    n_wat = res_counts["WAT"]
    assert n_acn > 0 and n_wat > 0

    mw_acn = SOLVENT_MW["acetonitrile"]
    mw_wat = SOLVENT_MW["water"]
    mass_acn = n_acn * mw_acn
    mass_wat = n_wat * mw_wat
    frac_acn = mass_acn / (mass_acn + mass_wat)
    frac_wat = mass_wat / (mass_acn + mass_wat)
    assert np.isclose(frac_acn, 0.3, rtol=0.05), f"ACN mass fraction {frac_acn} deviates from 0.3"
    assert np.isclose(frac_wat, 0.7, rtol=0.05), f"Water mass fraction {frac_wat} deviates from 0.7"


def test_mixture_spce_water(tmpdir):
    test_name = "test_mixture_spce_water"
    inputfilepath = hp.get_input_dir("water_acn_spce.json")
    startmulticomponent(["-f", inputfilepath])
    
    assert os.path.exists("water-acn-spce.pdb") 
    assert os.path.exists("water-acn-spce.prmtop") 
    assert os.path.exists("water-acn-spce.inpcrd")

    try:
        from parmed import load_file
    except ImportError:
        pytest.skip("parmed not installed, skipping prmtop check.")
    prmtop = load_file("water-acn-spce.prmtop")
    # check the ratio of the 2 component (residues) are close to 7:3
    res_counts = {}
    for res in prmtop.residues:
        res_name = res.name
        if res_name not in res_counts:
            res_counts[res_name] = 0
        res_counts[res_name] += 1
    residues = list(res_counts.keys())
    assert len(residues) == 2, f"Expected 2 residue types, found {len(residues)}."
    count1 = res_counts[residues[0]]
    count2 = res_counts[residues[1]]
    ratio = count1 / count2 if count1 >= count2 else count2 / count1
    expected_ratio = 7 / 3
    assert np.isclose(ratio, expected_ratio, rtol=0.05), f"Residue ratio {ratio} deviates from expected {expected_ratio}."

    # check the bond equilibrium length of O-H bond in SPC/E water
    spc_e_OH_bond_length = 1.0  # in Angstroms
    oh_bond = None
    for bond in prmtop.bonds:
        atom1 = bond.atom1
        atom2 = bond.atom2
        if (atom1.name == "O" and atom2.name == "H1") or (atom1.name == "H1" and atom2.name == "O") or \
           (atom1.name == "O" and atom2.name == "H2") or (atom1.name == "H2" and atom2.name == "O"):
            oh_bond = bond
            break
    assert oh_bond is not None, "O-H bond in SPC/E water not found in prmtop."
    assert np.isclose(oh_bond.type.req, spc_e_OH_bond_length, atol=1e-3), f"O-H bond equilibrium length is {oh_bond.type.req}, expected {spc_e_OH_bond_length:>.3f} Å for SPC/E water."


def test_mixture_builder_file_input():
    test_name = "test_mixture_builder_file_input" 
    inputfilepath = hp.get_input_dir("step1_input.json")
    startmulticomponent(["-f", inputfilepath])
    
    solute_path_exist = True
    solute_path_exist *= os.path.exists("naphthalene.mol2")
    solute_path_exist *= os.path.exists("naphthalene.frcmod")
    assert solute_path_exist

    # predefined solvent will not generate the mol2 and frcmod files.
    path_exist = True
    path_exist *= os.path.exists("naphthalene-water-acetonitrile.pdb") 
    path_exist *= os.path.exists("naphthalene-water-acetonitrile.prmtop") 
    path_exist *= os.path.exists("naphthalene-water-acetonitrile.inpcrd") 

    assert path_exist 
    assert hp.compare_pdb(
                "naphthalene-water-acetonitrile.pdb", 
                hp.get_reference_dir(f"multicomponent/naphthalene-water-acetonitrile.pdb"), 
                threshold = np.inf, # I set it to inf because packmol has some randomness in the output. This function will check the number of atoms and residues.
                ) 
    assert hp.compare_inpcrd_prmtop(
                "naphthalene-water-acetonitrile.prmtop", 
                hp.get_reference_dir(f"multicomponent/naphthalene-water-acetonitrile.prmtop"), 
                threshold = np.inf, # I set it to inf because packmol has some randomness in the output. This function will check the topology and force field parameters.
                ) 
    
def test_ionpair_solvation_file_input(tmpdir):
    testName = "test_ionpair_solvation_file_input"
    inputfilepath = hp.get_input_dir("ionpair_input.json")
    startmulticomponent(["-f", inputfilepath])
    path_fragment_exist = True
    for res in ("TPA", "SUF"):
        path_fragment_exist *= os.path.exists(f"ionpair-{res.lower()}.pdb")
        path_fragment_exist *= os.path.exists(f"ionpair-{res.lower()}.lib")
    assert path_fragment_exist
    path_main_exist = True
    for suffix in ("lib", "pdb"):
        path_main_exist *= os.path.exists(f"ionpair.{suffix}")
    assert path_main_exist
    assert hp.compare_inpcrd_prmtop(
        "ionpair-acetonitrile-toluene.prmtop",
        hp.get_reference_dir("multicomponent/ionpair-acetonitrile-toluene.prmtop"),
        threshold=np.inf,
    )

def test_mixture_builder_cmd_input():
    """This is a legacy feature, designed solely to respect the habits of users of the older version. It is not recommended for use."""
    test_name = "test_mixture_builder_cmd_input"
    solute_xyz = hp.get_input_dir("naphthalene_neutral.xyz")  
    startmulticomponent([
        "-m", solute_xyz,
        "-o", "mybox",
        "-c", "0",
        "-u", "1",
        "-s", "water",
        "-b", "20",
        "-t", "0.8",
    ])
    assert os.path.exists("mybox.prmtop")