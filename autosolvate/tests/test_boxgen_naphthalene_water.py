import autosolvate
import pytest
from . import helper_functions as hp
import os

def test_example_1(tmp_path):
    autosolvate.startboxgen(["-m", "inputs/naphthalene_neutral.xyz", "-o", "neutral"])
    compare_pdb, compare_inpcrd, compare_prmtop = compare_boxgen("neutral", os.path.join(hp.get_reference_dir(), "naphthalene_water/neutral"))
    assert compare_pdb
    assert compare_inpcrd
    assert compare_prmtop

    autosolvate.startboxgen(["-m", "inputs/naphthalene_radical.xyz", "-o", "radical"])
    compare_pdb, compare_inpcrd, compare_prmtop = compare_boxgen("radical", os.path.join(hp.get_reference_dir(), "naphthalene_water/radical"))
    assert compare_pdb
    assert compare_inpcrd
    assert compare_prmtop

#a general helper method to test pdb, prmtop and inpcrd files given output and reference filenames
def compare_boxgen(out, ref):
    with open(out + ".pdb", 'r') as f_in:
        pdb = f_in.read()
    with open(ref + ".pdb", 'r') as f_in:
        pdb_ref = f_in.read()
    compare_pdb = pdb == pdb_ref

    with open(out + ".inpcrd", 'r') as f_in:
        inpcrd = f_in.read()
    with open(ref + ".inpcrd", 'r') as f_in:
        inpcrd_ref = f_in.read()
    compare_inpcrd = inpcrd == inpcrd_ref

    with open(out + ".prmtop", 'r') as f_in:
        prmtop = f_in.readlines()[1:] #skip DATE = 10/08/22  11:32:06
    with open(ref + ".prmtop", 'r') as f_in:
        prmtop_ref = f_in.readlines()[1:]
    compare_prmtop = prmtop == prmtop_ref

    return compare_pdb, compare_inpcrd, compare_prmtop