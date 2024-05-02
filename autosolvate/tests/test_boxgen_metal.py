import os
import pytest

from autosolvate.boxgen_metal import *
from . import helper_functions as hp

def testboxgen_metal(tmpdir):
    os.system('cp inputs/QMresult/* ' + str(tmpdir))
    os.system('rm *.prmtop')
    os.system('rm *.inpcrd')
    os.system('rm *.cmd')
    os.system('rm *.log')
    startboxgen(["-m", "Fe_plus2", "-c", "0", "-s", "acetonitrile"])
    
    # 打印tleap.log的内容
    if os.path.exists('tleap.log'):
        with open('tleap.log', 'r') as file:
            print("Contents of tleap.log:")
            print(file.read())
    
    assert os.path.exists('tleap.log'), "tleap.log file does not exist."
    assert os.path.exists('Fe_plus2_solvated.prmtop'), "Fe_plus2_solvated.prmtop file does not exist."
    assert os.path.exists('Fe_plus2_solvated.inpcrd'), "Fe_plus2_solvated.inpcrd file does not exist."
