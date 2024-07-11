#--------------------------------------------------------------------------------------------------#
# test multicomponent amber parameter generation
# author: Fangning Ren (2022-12-17)
# path: autosolvate/tests/test_multicomponent.py 
#--------------------------------------------------------------------------------------------------#
import os
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

def test_ionpair_solvation_custom_solvent(tmpdir):
    testName = "test_custom_ionpair_solvation"
    solutexyz = hp.get_input_dir("ionpair.pdb")
    name = "ionpair"
    inst = MulticomponentSolventBoxBuilder(
        xyzfile = solutexyz,
        slu_charge={"SUF":-2, "TPA":1},
        solvent="acetonitrile",
        cube_size=50,
        charge_method="bcc",
        folder = os.getcwd(),
        custom_ionpair = True
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

def test_mixture_builder():
    test_name = "test_mixture_builder" 
    builder = MixtureBuilder(folder=os.getcwd()) 
    builder.add_solute(hp.get_input_dir("naphthalene_neutral.xyz"), name="naphthalene", residue_name="NAP", charge=0, spinmult=1, number=1)
   
    # predefined solvents do not need to specify the xyz file.
    builder.add_solvent(name="acetonitrile", residue_name="C3N", charge=0, spinmult=1, number=200)
    builder.add_solvent(name="water", residue_name="WAT", charge=0, spinmult=1, number=600)
    builder.build() 

    solute_path_exist = True
    solute_path_exist *= os.path.exists("naphthalene.mol2")
    solute_path_exist *= os.path.exists("naphthalene.frcmod")
    assert solute_path_exist

    # predefined solvent will not generate the mol2 and frcmod files.
    path_exist = True
    path_exist *= os.path.exists("naphthalene-acetonitrile-water.pdb") 
    path_exist *= os.path.exists("naphthalene-acetonitrile-water.prmtop") 
    path_exist *= os.path.exists("naphthalene-acetonitrile-water.inpcrd") 

    assert path_exist 
    assert hp.compare_pdb(
                "naphthalene-acetonitrile-water.pdb", 
                hp.get_reference_dir(f"multicomponent/naphthalene-acetonitrile-water.pdb"), 
                threshold = np.inf, # I set it to inf because packmol has some randomness in the output. This function will check the number of atoms and residues.
                ) 
    assert hp.compare_inpcrd_prmtop(
                "naphthalene-acetonitrile-water.prmtop", 
                hp.get_reference_dir(f"multicomponent/naphthalene-acetonitrile-water.prmtop"), 
                threshold = np.inf, # I set it to inf because packmol has some randomness in the output. This function will check the topology and force field parameters.
                ) 

def test_mixture_builder_file_input():
    test_name = "test_mixture_builder_file_input" 
    inputfilepath = hp.get_input_dir("mixturebuilder_input1.json")
    startmulticomponent(["-f", inputfilepath])
    
    solute_path_exist = True
    solute_path_exist *= os.path.exists("naphthalene.mol2")
    solute_path_exist *= os.path.exists("naphthalene.frcmod")
    assert solute_path_exist

    # predefined solvent will not generate the mol2 and frcmod files.
    path_exist = True
    path_exist *= os.path.exists("naphthalene-acetonitrile-water.pdb") 
    path_exist *= os.path.exists("naphthalene-acetonitrile-water.prmtop") 
    path_exist *= os.path.exists("naphthalene-acetonitrile-water.inpcrd") 

    assert path_exist 
    assert hp.compare_pdb(
                "naphthalene-acetonitrile-water.pdb", 
                hp.get_reference_dir(f"multicomponent/naphthalene-acetonitrile-water.pdb"), 
                threshold = np.inf, # I set it to inf because packmol has some randomness in the output. This function will check the number of atoms and residues.
                ) 
    assert hp.compare_inpcrd_prmtop(
                "naphthalene-acetonitrile-water.prmtop", 
                hp.get_reference_dir(f"multicomponent/naphthalene-acetonitrile-water.prmtop"), 
                threshold = np.inf, # I set it to inf because packmol has some randomness in the output. This function will check the topology and force field parameters.
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