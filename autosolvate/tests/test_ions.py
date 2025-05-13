import os
import pytest
import numpy as np

import autosolvate
from . import helper_functions as hp

def compare_boxgen(out, ref):
    with open(out + ".pdb", 'r') as f_in:
        pdb = f_in.read()
    with open(ref + ".pdb", 'r') as f_in:
        pdb_ref = f_in.read()
    compare_pdb = pdb == pdb_ref

    with open(out + ".inpcrd", 'r') as f_in:
        inpcrd = f_in.read()
    with open(ref + ".inpcrd", 'r') as f_in:
        inpcrd_ref = f_in.read()
    compare_inpcrd = inpcrd == inpcrd_ref

    with open(out + ".prmtop", 'r') as f_in:
        prmtop = f_in.readlines()[1:] #skip DATE = 10/08/22  11:32:06
    with open(ref + ".prmtop", 'r') as f_in:
        prmtop_ref = f_in.readlines()[1:]
    compare_prmtop = prmtop == prmtop_ref

    return compare_pdb, compare_inpcrd, compare_prmtop

def test_advanced_example_1(tmpdir):
    solventboxname = "water"
    os.system('cp inputs/K.xyz ' + str(tmpdir))
    
    
    
    autosolvate.startboxgen([
        "-m", "K.xyz",
        "-s", solventboxname,
        "-c", "1",
        "-o",'K_solvated'
        ])

    for p in tmpdir.listdir():
        print("  ", p.basename)

    pass_exist = True
    for suffix in ["K_solvated.pdb", "K_solvated.prmtop", "K_solvated.inpcrd"]:
        pass_exist *= os.path.exists(suffix)
    assert pass_exist
    
    
    compare_pdb, compare_inpcrd, compare_prmtop = compare_boxgen(os.path.join(tmpdir, "K_solvated"), os.path.join(hp.get_reference_dir(), "ions/K/K_solvated"))
    assert compare_pdb
    assert compare_inpcrd
    assert compare_prmtop
    
    os.system('cp inputs/S.xyz ' + str(tmpdir))
    
    solventboxname = "acetonitrile"
    autosolvate.startboxgen([
        "-m", "S.xyz",
        "-s", solventboxname,
        "-c", "-2",
        "-o",'S_solvated'
        ])

    for p in tmpdir.listdir():
        print("  ", p.basename)

    pass_exist = True
    for suffix in ["S_solvated.pdb", "S_solvated.prmtop", "S_solvated.inpcrd"]:
        pass_exist *= os.path.exists(suffix)
    assert pass_exist
    
    
    compare_pdb, compare_inpcrd, compare_prmtop = compare_boxgen(os.path.join(tmpdir, "S_solvated"), os.path.join(hp.get_reference_dir(), "ions/S/S_solvated"))
    assert compare_pdb
    assert compare_inpcrd
    assert compare_prmtop
    
    os.system('cp inputs/dmso.frcmod ' + str(tmpdir))
    os.system('cp inputs/dmso.off ' + str(tmpdir))
    os.system('cp inputs/Fe.xyz ' + str(tmpdir))
    
    solventboxname = "d"
    autosolvate.startboxgen([
        "-m", "Fe.xyz",
        "-s", solventboxname,
        "-c", "2",
        "-o",'Fe_dmso_solvated',
        "-l", "dmso.off",
        "-p", "dmso.frcmod"
        ])

    for p in tmpdir.listdir():
        print("  ", p.basename)

    pass_exist = True
    for suffix in ["Fe_dmso_solvated.pdb", "Fe_dmso_solvated.prmtop", "Fe_dmso_solvated.inpcrd"]:
        pass_exist *= os.path.exists(suffix)
    assert pass_exist

    compare_pdb, compare_inpcrd, compare_prmtop = compare_boxgen(os.path.join(tmpdir, "Fe_dmso_solvated"), os.path.join(hp.get_reference_dir(), "ions/Fe/Fe_dmso_solvated"))
    assert compare_pdb
    assert compare_inpcrd
    assert compare_prmtop