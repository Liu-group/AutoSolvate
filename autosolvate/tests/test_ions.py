import os
import pytest
import numpy as np

import autosolvate
from . import helper_functions as hp

def test_advanced_example_1(tmpdir):
    testName = "test_ions"
    solventboxname = "water"
    
    autosolvate.startboxgen([
        "-m", hp.get_input_dir("K.xyz"),
        "-s", solventboxname,
        "-c", "1"
        ])
    out = solventboxname + "_solvated"  
    ref = hp.get_reference_dir(solventboxname + "_solvated")

    pass_exist = True
    for suffix in [".pdb", ".prmtop", ".inpcrd"]:
        pass_exist *= os.path.exists(os.path.splitext(out)[0] + suffix)
    assert pass_exist

    pass_geom = hp.compare_pdb(out, ref)
    pass_amberinput = hp.compare_inpcrd_prmtop(out, ref)
    assert pass_geom
    assert pass_amberinput