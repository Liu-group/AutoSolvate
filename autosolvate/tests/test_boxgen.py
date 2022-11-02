import autosolvate
import pytest
import os
from pkg_resources import resource_filename, Requirement
from . import helper_functions as hp

def compare_boxgen(out, ref):
    with open(out + ".pdb", 'r') as f_in:
        pdb = f_in.read()
    with open(ref + ".pdb", 'r') as f_in:
        pdb_ref = f_in.read()

    with open(out + ".inpcrd", 'r') as f_in:
        inpcrd = f_in.read()
    with open(ref + ".inpcrd", 'r') as f_in:
        inpcrd_ref = f_in.read()

    with open(out + ".prmtop", 'r') as f_in:
        prmtop = f_in.readlines()[1:] #skip DATE
    with open(ref + ".prmtop", 'r') as f_in:
        prmtop_ref = f_in.readlines()[1:]

    return pdb == pdb_ref, inpcrd == inpcrd_ref, prmtop == prmtop_ref

def test_boxgen(tmpdir):
    autosolvate.startboxgen(["-m", hp.get_input_dir("naphthalene_neutral.xyz"),
                             "-s", "acetonitrile",
                             "-c", 0,
                             "-u", 1,
                             "-g", "bcc",
                             "-o", "nap_neutral_MeCN"])
    compare_exist = True
    for suffix in [".pdb", ".prmtop", ".inpcrd"]:
        compare_exist *= os.path.exists("nap_neutral_MeCN" + suffix)
    assert compare_exist 
    assert compare_boxgen("nap_neutral_MeCN",hp.get_reference_dir("nap_neutral_MeCN"))
