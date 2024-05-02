import os
import pytest
import numpy as np
import os

from autosolvate.boxgen_metal import *
from . import helper_functions as hp

def testboxgen_metal(tmpdir):
    os.system('cp inputs/QMresult/*.mol2 ' + tmpdir)
    os.system('cp inputs/QMresult/*.pdb ' + tmpdir)
    os.system('cp inputs/QMresult/*.frcmod ' + tmpdir)
    os.system('cp inputs/QMresult/*.info ' + tmpdir)
    os.system('cp inputs/QMresult/*.in ' + tmpdir)
    startboxgen(["-m", "Fe_plus2", "-c", "2","-u","1","-s","acetonitrile"])
    assert os.path.exists('Fe_plus2_solvated.prmtop')
    assert os.path.exists('Fe_plus2_solvated.inpcrd')