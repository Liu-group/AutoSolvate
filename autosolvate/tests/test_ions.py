import os
import pytest
import numpy as np

import autosolvate
from . import helper_functions as hp

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