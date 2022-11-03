import pytest 
from pathlib import Path
import os
from . import helper_functions as hp
import autosolvate

def test_temp_dir(tmpdir):
    """A function for testing the helperfunctions"""
    f = open(os.path.join(os.path.dirname(__file__), "log.log"), "w")

    tmpdirstr = str(tmpdir)
    f.write("\nTESTSETDIR:\t" + os.path.dirname(__file__))
    f.write("\nTMPDIR:\t" + tmpdirstr)
    f.write("\nCURRENTDIR:\t" + os.getcwd())
    assert tmpdirstr != os.path.dirname(__file__)
    assert tmpdirstr not in __file__
    assert tmpdirstr == os.getcwd()

    hp.get_input_dir()
    hp.get_reference_dir()
    f.write("\nINPUTDIR:\t" + hp.get_input_dir())
    f.write("\nTEMPINPUTDIR:\t" + os.path.dirname(hp.get_input_dir("water_solvated.prmtop")))
    f.write("\nREFERENCEDIR:\t" + hp.get_reference_dir())

    f.write("\nFILEININPUTS:\t" + str(sorted(os.listdir(hp.get_input_dir()))))
    f.write("\nTEMPINPUTS:\t" + str(sorted(os.listdir(os.path.join(os.path.dirname(__file__), "inputs")))))
    f.write("\nFILEINREFERENCE:\t" + str(os.listdir(hp.get_reference_dir())))
    assert os.path.exists(hp.get_input_dir("water_solvated.prmtop"))
    assert os.path.exists(hp.get_reference_dir("d_solvated.pdb"))
    assert os.path.exists(os.path.join(tmpdirstr, "inputs"))
    assert sorted(os.listdir(hp.get_input_dir())) != []
    assert sorted(os.listdir(hp.get_input_dir())) == sorted(os.listdir(os.path.join(os.path.dirname(__file__), "inputs")))
    f.close()