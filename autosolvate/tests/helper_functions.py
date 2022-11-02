import os
import random
import numpy as np
from pkg_resources import resource_filename, Requirement

# MDTraj is imported to help us read the amber input file.
# Autosolvate cannot run without mdtraj
import mdtraj as md

#a general helper method to test pdb, prmtop and inpcrd files given output and reference filenames

def initialize_random():
    """Set the random seed to zero to keep identical"""
    random.seed(0)
    np.random.seed(0)
initialize_random() # this function is instantly executed. 

def get_input_dir(name = ""):
    """will return the input directory if name is empty"""
    # use this directory if we move the folder "tests" to the Autosolvate-main directory with setup.py . 
    # return resource_filename(Requirement.parse("autosolvate"), "/tests/inputs/" + name)
    return resource_filename(Requirement.parse("autosolvate"), "autosolvate/tests/inputs/" + name)

def get_reference_dir(name = ""):
    """will return the reference directory if name is empty"""
    # use this directory if we move the folder "tests" to the Autosolvate-main directory with setup.py . 
    # return resource_filename(Requirement.parse("autosolvate"), "/tests/refs/" + name)
    return resource_filename(Requirement.parse("autosolvate"), "autosolvate/tests/refs/" + name)

def compare_pdb(out, ref, threshold = 1.0e-3):
    """Compare the two pdb file. A default threshold 1.0e-3 is set for fuzzy compare because the accuracy for the coordinates of pdb file is only 1.0e-3"""
    out, ref = os.path.splitext(out)[0] + ".pdb", os.path.splitext(ref)[0] + ".pdb"
    mol_out = md.load(out)
    mol_ref = md.load(ref)
    if mol_out.n_atoms != mol_ref.n_atoms:
        return False
    xyz_out, xyz_ref = mol_out.xyz[0], mol_ref.xyz[0]
    max_dist = np.max(np.sum((xyz_out - xyz_ref)**2, axis = 1)**0.5)
    return max_dist <= threshold

def compare_inpcrd_prmtop(out, ref, threshold = 1.0e-6):
    """Compare the atomic RMSD for the inpcrd file for md simulation. Simultaneously compare the inpcrd file with coordinates and the prmtop file with charge and topology. They should have the exact identical topology."""
    out_prmtop, ref_prmtop = os.path.splitext(out)[0] + ".prmtop", os.path.splitext(ref)[0] + ".prmtop"
    out_inpcrd, ref_inpcrd = os.path.splitext(out)[0] + ".inpcrd", os.path.splitext(ref)[0] + ".inpcrd"
    mol_out = md.load(out_inpcrd, top = out_prmtop)
    mol_ref = md.load(ref_inpcrd, top = ref_prmtop)
    if mol_out.n_atoms != mol_ref.n_atoms:
        return False
    xyz_out, xyz_ref = mol_out.xyz[0], mol_ref.xyz[0]
    rmsd = (np.sum((xyz_out - xyz_ref)**2) / mol_out.n_atoms)**0.5
    if rmsd >= threshold:
        return False
    bonds_out, bonds_ref = list(mol_out.top.bonds), list(mol_ref.top.bonds)
    if len(bonds_out) != len(bonds_ref):
        return False
    idxs = np.random.randint(0, len(bonds_out), size = max(len(bonds_out), 50))
    for i in idxs:
        if bonds_out[i].atom1 != bonds_ref[i].atom1:
            return False
        if bonds_out[i].atom2 != bonds_ref[i].atom2:
            return False
    return True
    
def compare_boxgen(out, ref):
    compareExist = True
    for suffix in [".pdb", ".prmtop", ".inpcrd"]:
        compare_exist *= os.path.exists(os.path.splitext(out)[0] + suffix)
    comparePdb = compare_pdb(out, ref)
    compareInpcrd, comparePrmtop = compare_pdb(out, ref)
    return comparePdb, compareInpcrd, comparePrmtop