#--------------------------------------------------------------------------------------------------#
# Helper functions. A set of functions that can get the file path, read and compare the xyz and pdb files.
# author: Fangning Ren (2022-11-03) 
# path: autosolvate/tests/helper_functions.py
#--------------------------------------------------------------------------------------------------#
from collections import Counter
from pkg_resources import resource_filename, Requirement, DistributionNotFound
import os
import random

import numpy as np
import pytest
import parmed as pmd


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
    """
    # use this directory if we move the folder "tests" to the Autosolvate-main directory with setup.py . 
    path1 = os.path.join(os.getcwd(), "inputs")
    path2 = os.path.join(os.path.dirname(__file__), "inputs")
    try:
        path3 = resource_filename(Requirement.parse("autosolvate"), "autosolvate/tests/inputs")
    except DistributionNotFound:
        path3 = path2
    for ppath in [path1, path2, path3]:
        if not os.path.exists(ppath):
            continue
        if os.path.exists(os.path.join(ppath, name)):
            return os.path.join(ppath, name)
        flist = os.listdir(ppath)
        flist = [fname for fname in flist if fname.find(name) != -1]
        if len(flist) >= 1:
            return os.path.join(ppath, name)
    raise FileNotFoundError(os.path.join(path1, name) + "\t" + os.path.join(path2, name) + "\t not exist.")

def get_reference_dir(name = ""):
    """will return the reference directory if name is empty"""
    path2 = os.path.join(os.path.dirname(__file__), "refs")
    try:
        path3 = resource_filename(Requirement.parse("autosolvate"), "autosolvate/tests/refs")
    except DistributionNotFound:
        path3 = path2
    for ppath in [path2, path3]:
        if not os.path.exists(ppath):
            continue
        if os.path.exists(os.path.join(ppath, name)):
            return os.path.join(ppath, name)
        flist = os.listdir(ppath)
        flist = [fname for fname in flist if fname.find(name) != -1]
        if len(flist) >= 1:
            return os.path.join(ppath, name)
    raise FileNotFoundError(os.path.join(path2, name) + "\t not exist.")

def get_temporary_dir(tmpdir, name = ""):
    """will return the temporary directory if name is empty. Can be replaced by tmpdir.dirname + '/' + name"""
    return str(tmpdir / name)

# def compare_pdb(out, ref, threshold = 1.0e-3):
#     """Compare the two pdb file. A default threshold 1.0e-3 is set for fuzzy compare because the accuracy for the coordinates of pdb file is only 1.0e-3"""
#     out, ref = os.path.splitext(out)[0] + ".pdb", os.path.splitext(ref)[0] + ".pdb"
#     mol_out = md.load_pdb(out)
#     mol_ref = md.load_pdb(ref)
#     if mol_out.n_atoms != mol_ref.n_atoms:
#         return False
#     top_out, top_ref = mol_out.top, mol_ref.top
#     top_out:md.Topology
#     for i in range(mol_out.n_atoms):
#         if top_out.atom(i).name != top_ref.atom(i).name:
#             return False
#         if top_out.atom(i).residue.name != top_ref.atom(i).residue.name:
#             return False
#     xyz_out, xyz_ref = mol_out.xyz[0], mol_ref.xyz[0]
#     max_dist = np.max(np.sum((xyz_out - xyz_ref)**2, axis = 1)**0.5)
#     return max_dist <= threshold

def compare_pdb(out, ref, threshold = 1.0e-3):
    """Compare the two pdb file. A default threshold 1.0e-3 is set for fuzzy compare because the accuracy for the coordinates of pdb file is only 1.0e-3"""
    out, ref = os.path.splitext(out)[0] + ".pdb", os.path.splitext(ref)[0] + ".pdb"
    mol_out = pmd.load_file(out)
    mol_ref = pmd.load_file(ref)
    if len(mol_out.atoms) != len(mol_ref.atoms):
        return False
    for i in range(len(mol_out.atoms)):
        if mol_out.atoms[i].name != mol_ref.atoms[i].name:
            return False
        if mol_out.atoms[i].residue.name != mol_ref.atoms[i].residue.name:
            return False
    xyz_out, xyz_ref = mol_out.coordinates, mol_ref.coordinates
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

# def compare_inpcrd_prmtop(out, ref, threshold = 1.0e-6):
#     """Compare the atomic RMSD for the inpcrd file for md simulation. Simultaneously compare the inpcrd file with coordinates and the prmtop file with charge and topology. They should have the exact identical topology."""
#     out_prmtop, ref_prmtop = os.path.splitext(out)[0] + ".prmtop", os.path.splitext(ref)[0] + ".prmtop"
#     out_inpcrd, ref_inpcrd = os.path.splitext(out)[0] + ".inpcrd", os.path.splitext(ref)[0] + ".inpcrd"
#     mol_out = md.load(out_inpcrd, top = out_prmtop)
#     mol_ref = md.load(ref_inpcrd, top = ref_prmtop)
#     if mol_out.n_atoms != mol_ref.n_atoms:
#         return False
#     xyz_out, xyz_ref = mol_out.xyz[0], mol_ref.xyz[0]
#     rmsd = (np.sum((xyz_out - xyz_ref)**2) / mol_out.n_atoms)**0.5
#     if rmsd >= threshold:
#         return False
#     bonds_out, bonds_ref = list(mol_out.top.bonds), list(mol_ref.top.bonds)
#     if len(bonds_out) != len(bonds_ref):
#         return False
#     idxs = np.random.randint(0, len(bonds_out), size = max(len(bonds_out), 50))
#     for i in idxs:
#         if bonds_out[i].atom1 != bonds_ref[i].atom1:
#             return False
#         if bonds_out[i].atom2 != bonds_ref[i].atom2:
#             return False
#     return True

def compare_inpcrd_prmtop(out, ref, threshold = 1.0e-6):
    """Compare the atomic RMSD for the inpcrd file for md simulation. Simultaneously compare the inpcrd file with coordinates and the prmtop file with charge and topology. They should have the exact identical topology."""
    out_prmtop, ref_prmtop = os.path.splitext(out)[0] + ".prmtop", os.path.splitext(ref)[0] + ".prmtop"
    out_inpcrd, ref_inpcrd = os.path.splitext(out)[0] + ".inpcrd", os.path.splitext(ref)[0] + ".inpcrd"
    mol_out = pmd.load_file(out_prmtop, xyz = out_inpcrd)
    mol_ref = pmd.load_file(ref_prmtop, xyz = ref_inpcrd)
    if len(mol_out.atoms) != len(mol_ref.atoms):
        print(len(mol_out.atoms), len(mol_ref.atoms))
        return False
    xyz_out, xyz_ref = mol_out.coordinates, mol_ref.coordinates
    rmsd = (np.sum((xyz_out - xyz_ref)**2) / len(mol_out.atoms))**0.5
    if rmsd >= threshold:
        print(rmsd)
        return False
    bonds_out, bonds_ref = mol_out.bonds, mol_ref.bonds
    if len(bonds_out) != len(bonds_ref):
        print(len(bonds_out), len(bonds_ref))
        return False
    idxs = np.random.randint(0, len(bonds_out), size = max(len(bonds_out), 50))
    for i in idxs:
        if (bonds_out[i].atom1.name != bonds_ref[i].atom1.name or bonds_out[i].atom2.name != bonds_ref[i].atom2.name) and \
           (bonds_out[i].atom1.name != bonds_ref[i].atom2.name or bonds_out[i].atom2.name != bonds_ref[i].atom1.name): 
            print(bonds_out[i].atom1.name, bonds_ref[i].atom1.name, bonds_out[i].atom2.name, bonds_ref[i].atom2.name)
            return False
    return True
    
def compare_boxgen(out, ref):
    compareExist = True
    for suffix in [".pdb", ".prmtop", ".inpcrd"]:
        compareExist *= os.path.exists(os.path.splitext(out)[0] + suffix)
    comparePdb = compare_pdb(out, ref)
    compareInpcrd, comparePrmtop = compare_inpcrd_prmtop(out, ref)
    return comparePdb, compareInpcrd, comparePrmtop


def check_coordinate_number(prmtoppath, metal_residue_name, coordinate_number = 6):
    prmtop = pmd.amber.AmberParm(prmtoppath)
    metal_legand_bonds = []
    for bond in prmtop.bonds:
        a1, a2 = bond.atom1, bond.atom2
        if a1.residue.name == metal_residue_name and a2.residue.name != metal_residue_name:
            metal_legand_bonds.append(bond)
        elif a2.residue.name == metal_residue_name and a1.residue.name != metal_residue_name:
            metal_legand_bonds.append(bond)
    return len(metal_legand_bonds) == coordinate_number