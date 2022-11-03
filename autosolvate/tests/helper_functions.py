#--------------------------------------------------------------------------------------------------#
# Helper functions. A set of functions that can get the file path, read and compare the xyz and pdb files.
# author: Fangning Ren (2022-11-03) 
# path: autosolvate/tests/helper_functions.py
#--------------------------------------------------------------------------------------------------#
import os
import random
import numpy as np
from pkg_resources import resource_filename, Requirement, DistributionNotFound
from pathlib import *
import pytest
import autosolvate

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
    """
    will return the input directory if name is empty
    This function will first check the temporary directory. If it do not exist, find it at the input directory.
    This function can accept filename without extension. It will just add the absolute path before the file name.
    """
    # use this directory if we move the folder "tests" to the Autosolvate-main directory with setup.py . 
    if os.path.exists(os.path.join(os.getcwd(), "inputs")):
        return os.path.join(os.getcwd(), "inputs", name)
    if os.path.exists(os.path.join(os.path.dirname(__file__), "inputs")):
        return os.path.join(os.path.dirname(__file__), "inputs", name)
    try:
        return resource_filename(Requirement.parse("autosolvate"), "autosolvate/tests/inputs/" + name)
    except DistributionNotFound:
        pass
    raise FileNotFoundError(os.path.join(os.path.dirname(__file__), "inputs", name) + "\t" + os.path.join(os.getcwd(), "inputs", name) + "\t not exist.")

def get_reference_dir(name = ""):
    """will return the reference directory if name is empty"""
    # use this directory if we move the folder "tests" to the Autosolvate-main directory with setup.py . 
    if os.path.exists(os.path.join(os.path.dirname(__file__), "refs")):
        return os.path.join(os.path.dirname(__file__), "refs", name)
    try:
        return resource_filename(Requirement.parse("autosolvate"), "autosolvate/tests/refs/" + name)
    except DistributionNotFound:
        pass
    raise FileNotFoundError(os.path.join(os.path.dirname(__file__), "refs", name)+"\t not exist.")

def get_temporary_dir(tmpdir:Path, name = ""):
    """will return the temporary directory if name is empty. Can be replaced by tmpdir.dirname + '/' + name"""
    return str(tmpdir / name)

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

def read_xyz_multiple(fname):
    """simple function for reading a xyz file contains multiple conformations"""
    f = open(fname, "r")
    lines = f.readlines()
    f.close()
    curlineidx = 0
    elements, datas = [], []
    while curlineidx < len(lines) and lines[curlineidx].strip():
        natom = int(lines[curlineidx])
        datastr = lines[curlineidx+2:curlineidx+2+natom]
        data = np.empty((natom, 3), dtype = float)
        element = []
        for i, line in enumerate(datastr):
            temp = line.split()
            if len(temp) != 4:
                continue
            a, x, y, z = temp
            element.append(a)
            data[i][0], data[i][1], data[i][2] = float(x), float(y), float(z)
        elements.append(element)
        datas.append(data)
        curlineidx = curlineidx+2+natom
    return elements, datas

def compare_xyz(out, ref, threshold = 1.0e-6):
    """Compare the two xyz file. A default threshold 1.0e-3 is set for fuzzy compare"""
    out, ref = os.path.splitext(out)[0] + ".xyz", os.path.splitext(ref)[0] + ".xyz"
    elems_out, xyz_out = read_xyz_multiple(out)
    elems_ref, xyz_ref = read_xyz_multiple(ref)
    elems_out, xyz_out = elems_out[0], xyz_out[0]
    elems_ref, xyz_ref = elems_ref[0], xyz_ref[0]
    if elems_out != elems_ref:
        return False
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