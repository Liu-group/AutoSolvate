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
    pass_exist = True
    for res in ("UAA", "UAB", "UAC", "UAD"):
        pass_exist *= os.path.exists(f"{inpfname}-{res.lower()}.pdb")
        pass_exist *= os.path.exists(f"{inpfname}-{res.lower()}.lib")
        pass_exist *= os.path.exists(f"{inpfname}-{res.lower()}.frcmod")
    assert pass_exist
    assert hp.compare_pdb(f"{inpfname}.pdb", hp.get_reference_dir(f"multicomponent/{inpfname}-processed.pdb"))

def test_ionpair_solvation(tmpdir):
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
    pass_fragment_exist = True
    for res in ("TPA", "SUF"):
        pass_fragment_exist *= os.path.exists(f"{name}-{res.lower()}.pdb")
        pass_fragment_exist *= os.path.exists(f"{name}-{res.lower()}.lib")
    assert pass_fragment_exist
    pass_main_exist = True
    for suffix in ("lib", "pdb"):
        pass_main_exist *= os.path.exists(f"{name}.{suffix}")
    assert pass_main_exist

    assert hp.compare_pdb(f"water_solvated.pdb", hp.get_reference_dir(f"multicomponent/water_solvated.pdb"))
    assert hp.compare_inpcrd_prmtop(f"water_solvated.prmtop", hp.get_reference_dir(f"multicomponent/water_solvated.prmtop"))


