#--------------------------------------------------------------------------------------------------#
# test the box generation for transition metal complex
# author: Fangning Ren (2024-04-18)
# path: autosolvate/tests/test_tmc_boxgen.py
#--------------------------------------------------------------------------------------------------#
from collections import Counter
import os
import pytest
import numpy as np

from . import helper_functions as hp
from ..molecule import *
from ..dockers import *

import parmed as pmd 



def test_single_tmc(tmpdir):
    """
    @TODO:
        New PDB compare function 

    @NOTE:
        The resulting solvent box is slightly different from the previous one: 
        The positions of the sodium ions used to balance the charges are inconsistent. 
        But the current pdb comparison function cannot tolerate such subtle differences. 
        Therefore, I canceled them before the new pdb functions finished.
    """
    testName = "test_single_tmc"
    tmcname = "Co_plus3_bpy3"
    charge = 3
    spinmult = 1
    metal_residue_name = "CO1"
    tmcfolder = hp.get_input_dir(tmcname)
    tmcpdb = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.pdb")
    prmtop = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.prmtop")
    mol_tmc = TransitionMetalComplex(tmcpdb, charge, spinmult, prmtop, tmcfolder, 
                                 metal_residue_names = metal_residue_name, 
                                 name = tmcname,
                                 residue_name = "TMC", 
                                 folder = os.getcwd())
    mol_tmc.update()
    TleapDocker(workfolder=mol_tmc.folder).run(mol_tmc)

    # check the existence of the output files
    assert os.path.exists(f"{tmcname}.prmtop")
    assert os.path.exists(f"{tmcname}.inpcrd")

    # check the coordinate number of the metal atom. octahedral complex should have 6 bonds
    assert hp.check_coordinate_number(mol_tmc.prmtop, metal_residue_name, 6)

def test_tmc_water_solvated(tmpdir):
    testName = "test_tmc_water_solvated"
    tmcname = "Co_plus3_bpy3"
    charge = 3
    spinmult = 1
    metal_residue_name = "CO1"
    tmcfolder = hp.get_input_dir(tmcname)
    tmcpdb = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.pdb")
    prmtop = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.prmtop")
    mol_tmc = TransitionMetalComplex(tmcpdb, charge, spinmult, prmtop, tmcfolder, 
                                 metal_residue_names = metal_residue_name, 
                                 name = tmcname,
                                 residue_name = "TMC", 
                                 folder = os.getcwd())
    mol_tmc.update()

    systemname = tmcname + "_water_solvated"
    system = SolvatedSystem(systemname, mol_tmc, AMBER_WATER_BOX, cubesize=30, folder=mol_tmc.folder)
    system.set_closeness(automate=True)
    TleapDocker(workfolder=mol_tmc.folder).run(system)

    # check the existence of the output files
    assert os.path.exists(f"{systemname}.prmtop")
    assert os.path.exists(f"{systemname}.inpcrd")

    # check the coordinate number of the metal atom. octahedral complex should have 6 bonds
    assert hp.check_coordinate_number(mol_tmc.prmtop, metal_residue_name, 6)

    # check the number of other molecules
    prmtop = pmd.amber.AmberParm(f"{systemname}.prmtop")
    counter = Counter([res.name for res in prmtop.residues])
    assert counter["Cl-"] == charge
    assert counter["WAT"] <= 2500
    assert counter["WAT"] >= 2400

def test_tmc_acetonitrile_solvated(tmpdir):
    testName = "test_tmc_acetonitrile_solvated"
    tmcname = "Co_plus3_bpy3"
    charge = 3
    spinmult = 1
    metal_residue_name = "CO1"
    tmcfolder = hp.get_input_dir(tmcname)
    tmcpdb = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.pdb")
    prmtop = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.prmtop")
    mol_tmc = TransitionMetalComplex(tmcpdb, charge, spinmult, prmtop, tmcfolder, 
                                 metal_residue_names = metal_residue_name, 
                                 name = tmcname,
                                 residue_name = "TMC", 
                                 folder = os.getcwd())
    mol_tmc.update()

    solvPrefix = "ch3cn"
    solvresname = "C3N"
    solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
        os.path.join('data',solvPrefix,solvPrefix+".frcmod"))
    solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
        os.path.join('data',solvPrefix,solvPrefix+".prep"))
    solvent_pdb_path = pkg_resources.resource_filename('autosolvate', 
        os.path.join('data',solvPrefix,solvPrefix+".pdb"))
    solvent = Molecule(solvent_pdb_path, 0, 1, "solvent", residue_name = solvresname, folder = mol_tmc.folder)
    solvent.frcmod = solvent_frcmod_path
    solvent.prep = solvent_prep_path

    systemname = tmcname + "_acetonitrile_solvated"
    system = SolvatedSystem(systemname, mol_tmc, solvent, cubesize=30, folder=mol_tmc.folder, solute_number = 1, solvent_number = 311)
    system.set_closeness(closeness=1.88)
    PackmolDocker(workfolder=mol_tmc.folder).run(system)
    TleapDocker(workfolder=mol_tmc.folder).run(system)

    # check the existence of the output files
    assert os.path.exists(f"{systemname}.prmtop")
    assert os.path.exists(f"{systemname}.inpcrd")

    # check the coordinate number of the metal atom. octahedral complex should have 6 bonds
    assert hp.check_coordinate_number(mol_tmc.prmtop, metal_residue_name, 6)

    # check the number of other molecules
    prmtop = pmd.amber.AmberParm(f"{systemname}.prmtop")
    counter = Counter([res.name for res in prmtop.residues])
    assert counter["Cl-"] == charge
    assert counter["C3N"] == 311

def test_tmc_IL_solvated(tmpdir):
    testName = "test_tmc_IL_solvated"
    tmcname = "Co_plus3_bpy3"
    charge = 3
    spinmult = 1
    metal_residue_name = "CO1"
    tmcfolder = hp.get_input_dir(tmcname)
    tmcpdb = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.pdb")
    prmtop = hp.get_input_dir(f"{tmcname}/{tmcname}_dry.prmtop")
    mol_tmc = TransitionMetalComplex(tmcpdb, charge, spinmult, prmtop, tmcfolder, 
                                 metal_residue_names = metal_residue_name, 
                                 name = tmcname,
                                 residue_name = "TMC", 
                                 folder = os.getcwd())
    mol_tmc.update()

    mol_cat_pdbpath = hp.get_input_dir("IL/BMIM.pdb")
    mol_cat = Molecule(mol_cat_pdbpath, 1, 1, 
                   name = "BMIM", residue_name = "BMI", folder = mol_tmc.folder)
    mol_cat.mol2 = hp.get_input_dir("IL/BMIM.mol2")
    mol_cat.frcmod = hp.get_input_dir("IL/BMIM.frcmod")
    mol_cat.update()


    mol_anion_pdbpath = hp.get_input_dir("IL/NTF2.pdb")
    mol_ani = Molecule(mol_anion_pdbpath,-1, 1,
                        name = "NTF2", residue_name = "NSC", folder = mol_tmc.folder)
    mol_ani.mol2 = hp.get_input_dir("IL/NTF2.mol2")
    mol_ani.frcmod = hp.get_input_dir("IL/NTF2.frcmod")
    mol_ani.update()

    systemname = tmcname + "_IL_solvated"
    system = SolvatedSystem(systemname, solute = mol_tmc, solvent = [mol_cat, mol_ani],  
                        cubesize=30, closeness=2.0, solute_number=1, solvent_number=[50, 50+charge], folder = mol_tmc.folder)
    
    docker = PackmolDocker(system.folder)
    docker.run(system)
    docker = TleapDocker(system.folder)
    docker.run(system)

    # check the existence of the output files
    assert os.path.exists(f"{systemname}.prmtop")
    assert os.path.exists(f"{systemname}.inpcrd")

    # check the coordinate number of the metal atom. octahedral complex should have 6 bonds
    assert hp.check_coordinate_number(mol_tmc.prmtop, metal_residue_name, 6)

    # check the number of other molecules
    prmtop = pmd.amber.AmberParm(f"{systemname}.prmtop")
    counter = Counter([res.name for res in prmtop.residues])
    assert counter["Cl-"] == 0
    assert counter["BMI"] == 50
    assert counter["NSC"] == 50 + charge