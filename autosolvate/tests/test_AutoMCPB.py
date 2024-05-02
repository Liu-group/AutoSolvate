import os
import pytest
import numpy as np
import os

from autosolvate.autoMCPB import *
from . import helper_functions as hp


def testAutoMCPB_step1(tmpdir):
    os.system('cp inputs/QMresult/Fe_plus2.xyz ' + str(tmpdir))
    startautoMCPB(["-n", "Fe_plus2", "-c", "2","-u","1","-s","1","-x","orca"])
    assert os.path.exists('LG0.mol2')
    assert os.path.exists('LG1.mol2')
    assert os.path.exists('LG2.mol2')
    assert os.path.exists('LG0.frcmod')
    assert os.path.exists('LG1.frcmod')
    assert os.path.exists('LG2.frcmod')
    assert os.path.exists('FE.mol2')
    assert os.path.exists('Fe_plus2_final.pdb')
    assert os.path.exists('Fe_plus2_MCPB.in')


def testAutoMCPB_steps(tmpdir):
    os.system('cp inputs/QMresult/* ' + str(tmpdir))
    os.system('rm Fe_plus2_mcpbpy.frcmod')
    os.system('rm resp2.chg')
    os.system('rm Fe_plus2_dry.prmtop')
    os.system('ls *')

    startautoMCPB(["-n", "Fe_plus2", "-c", "2","-u","1","-s","2","-x","orca"])
    assert os.path.exists('Fe_plus2_mcpbpy.frcmod')
   # assert os.path.exists('Fe_plus2_mcpbpy_pre.frcmod')

    startautoMCPB(["-n", "Fe_plus2", "-c", "2","-u","1","-s","3","-x","orca"])
    assert os.path.exists('resp2.chg')

    startautoMCPB(["-n", "Fe_plus2", "-c", "2","-u","1","-s","4","-x","orca"])
    assert os.path.exists('Fe_plus2_dry.prmtop')





