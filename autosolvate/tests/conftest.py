#--------------------------------------------------------------------------------------------------#
# Conftest. These functions will automatically run by pytest.
# author: Fangning Ren (2022-11-03) 
# path: autosolvate/tests/helper_functions.py
#--------------------------------------------------------------------------------------------------#
import os
import shutil
import pytest
import subprocess
from . import helper_functions as hp

@pytest.fixture(autouse=True, scope = "function")
def _change2tmpdir(request):
    """This function will run before each test. The directory will be changed into the temporary directory"""
    tmpdir = request.getfixturevalue("tmpdir")
    with tmpdir.as_cwd():
        yield

@pytest.fixture(autouse=True, scope = "function")
def _copy_inputs(request):
    """All input files are copied into a folder called "input" in temporary directory"""
    # It is also dangerous to put all temporary output and input files together at the temporary directory
    # Therefore create an inputs directory to store the input files.
    inputpath = os.path.join(os.path.dirname(__file__), "inputs")
    inputfiles = os.listdir(inputpath)
    tmpdir = request.getfixturevalue("tmpdir")
    tmpdirstr = str(tmpdir)
    os.makedirs(os.path.join(tmpdirstr, "inputs"), exist_ok=True)
    for inputfile in inputfiles:
        if os.path.isdir(os.path.join(inputpath, inputfile)):
            shutil.copytree(os.path.join(inputpath, inputfile), os.path.join(tmpdirstr, "inputs", inputfile))
        else:
            shutil.copy(os.path.join(inputpath, inputfile), os.path.join(tmpdirstr, "inputs"))
    
@pytest.fixture(autouse=True, scope = "session")
def _activate_autosolvate():
    return 
    # result = subprocess.run("conda list")
