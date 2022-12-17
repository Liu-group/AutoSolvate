#--------------------------------------------------------------------------------------------------#
# test multicomponent amber parameter generation
# author: Fangning Ren (2022-12-17)
# path: autosolvate/tests/test_multicomponent.py 
#--------------------------------------------------------------------------------------------------#
import os
import pytest
import numpy as np

import autosolvate
from autosolvate.multicomponent import MulticomponentParamsBuilder
from . import helper_functions as hp


def test_multicomponent(tmpdir):
    testName = "test_multicomponent"
    inpfname = "PAHs"
    builder = MulticomponentParamsBuilder(
        hp.get_input_dir(f"{inpfname}.pdb")
    )
    builder.buildAmberParamsForAll()
    pass_exist = True
    for res in ("UAA", "UAB", "UAC", "UAD"):
        pass_exist *= os.path.exists(f"{inpfname}-{res.lower()}.pdb")
        pass_exist *= os.path.exists(f"{inpfname}-{res.lower()}.lib")
    assert pass_exist
    assert os.path.exists(f"{inpfname}.pdb")
    assert os.path.exists(f"{inpfname}.lib")
    assert hp.compare_pdb(f"{inpfname}.pdb", hp.get_reference_dir(f"{inpfname}-processed.pdb"))

