import os
import subprocess
import shutil
import numpy as np
from copy import deepcopy
from typing import List, Tuple, Dict, Optional, Union, TextIO

#import mdtraj as md

from ..molecule import *
from .general_docker import GeneralDocker
from ..utils import runMM, read_qm_folder

"""
TODO: 
    1. Finish check function
    2. correct the logging behavior
"""

class ForceBalanceDocker(GeneralDocker):
    def __init__(self, workfolder:str="wlp"):
        super().__init__(executable="ForceBalance", workfolder=workfolder)
        
        self.forcefield_dir = os.path.join(self.workfolder, "forcefield")
        self.result_dir     = os.path.join(self.workfolder, "result")
        self.targets_dir    = os.path.join(self.workfolder, "targets")
        self.cluster_dirs   = []
        self.n_clusters     = 0
        
        os.makedirs(self.forcefield_dir, exist_ok=True)
        os.makedirs(self.result_dir, exist_ok=True)
        os.makedirs(self.targets_dir, exist_ok=True)

        self.target_frcmod  = ""
        self.target_mol2    = ""
        self.result_frcmod  = ""

        self.inpfile = os.path.join(self.workfolder, "inpfile.inp")
        self.outfile = os.path.join(self.workfolder, "outfile.out")

    def write_qdata(self, qdatafile:str, xyzs, grds, enes):
        fo = open(qdatafile, "w")
        nconf = len(xyzs)
        for i in range(nconf):
            xyz, grd, ene = xyzs[i], grds[i], enes[i]
            fo.write(f"JOB {i}\n")
            fo.write(f"COORDS " + " ".join(["{:>17.10e}".format(a) for a in xyz.ravel()]) + "\n")
            fo.write(f"ENERGY " + "{:>17.10e}".format(ene) + "\n")
            fo.write(f"FORCES " + " ".join(["{:>17.10e}".format(a) for a in grd.ravel()]) + "\n")
            fo.write("\n")
        fo.close()

    def write_leap_setup(self, mol:Molecule, cluster_dir:str):
        with open(os.path.join(cluster_dir, "setup.leap"), "w") as fo:
            fo.write("source leaprc.gaff            \n")
            fo.write("source leaprc.protein.ff14SB  \n")
            fo.write("source leaprc.water.tip3p     \n")
            fo.write("loadamberparams {}            \n".format(os.path.basename(mol.frcmod)))
            fo.write("{} = loadmol2 {}              \n".format(mol.residue_name, os.path.basename(mol.mol2)))
            fo.write("check {}                      \n".format(mol.residue_name))
            fo.write("savepdb {} {}.pdb             \n".format(mol.residue_name, mol.name))
            fo.write("saveamberparm {} {}.prmtop {}.inpcrd\n".format(mol.residue_name, mol.name, mol.name))
            fo.write("quit\n")

    def add_data(self, mol:Molecule, qmfolder:str, prefix:str = ""):
        prefix = prefix if prefix else mol.name
        xyzs, grds, enes = read_qm_folder(qmfolder, prefix)
        nconf = len(xyzs)
        if nconf == 0:
            print("No data found for", mol.name)
            return
        curclusterid = self.n_clusters
        cluster_dir = os.path.join(self.targets_dir, "cluster-{:02d}".format(curclusterid))
        os.makedirs(cluster_dir, exist_ok=True)
        qdatafile = os.path.join(cluster_dir, f"qdata.txt")
        self.write_qdata(qdatafile, xyzs, grds, enes)

        t0 = md.load_xyz(mol.xyz, top = mol.prmtop)
        cfms = []
        for i in range(nconf):
            t = deepcopy(t0)
            t.xyz[0] = xyzs[i] / 10.0
            cfms.append(t)
        traj = md.join(cfms)
        traj.save_mdcrd(os.path.join(cluster_dir, f"all.mdcrd"))
        traj.slice(0).save_pdb(os.path.join(cluster_dir, f"conf.pdb"))

        self.write_leap_setup(mol, cluster_dir)

        self.cluster_dirs.append(cluster_dir)
        self.n_clusters += 1

    def gen_forcefield(self, mol:Molecule,):
        shutil.copy(mol.mol2, self.forcefield_dir)
        fmod = FrcmodFile(mol.frcmod)
        frcmodpath = os.path.join(self.forcefield_dir, os.path.basename(mol.frcmod))
        fmod.write_with_prm_label(frcmodpath)

    def check_system(self, mol:Molecule):
        pass

    def predict_output(self, mol:Molecule):
        self.inpfile        = os.path.join(self.workfolder, "{}-wlpfit.inp" .format(mol.name))
        self.outfile        = os.path.join(self.workfolder, "{}-wlpfit.out" .format(mol.name))
        self.target_frcmod  = os.path.join(self.targets_dir,"{}.frcmod"     .format(mol.name))
        self.target_mol2    = os.path.join(self.targets_dir,"{}.mol2"       .format(mol.name))
        self.result_frcmod  = os.path.join(self.result_dir, "{}.frcmod"     .format(mol.name))

    def write_cluster(self, fo:TextIO, cluster_dir:str):
        cluster_name = os.path.relpath(cluster_dir, self.targets_dir)
        fo.write(f"$target            \n")
        fo.write(f"w_energy 1.0       \n")
        fo.write(f"w_force 1.0        \n")
        fo.write(f"type abinitio_amber\n")
        fo.write(f"name {cluster_name}\n")
        fo.write(f"$end\n\n")
        
    def generate_input(self, mol:Molecule):
        frcmod_basename = os.path.basename(mol.frcmod)
        mol2_basename = os.path.basename(mol.mol2)
        fo = open(self.inpfile, "w")
        fo.write("$options              \n")
        fo.write("jobtype newton        \n")
        fo.write("forcefield {} {}      \n".format(frcmod_basename, mol2_basename))
        fo.write("penalty_additive 0.001\n")
        fo.write("trust0 1.0            \n")
        fo.write("$end                \n\n")
        for cluster_dir in self.cluster_dirs:
            self.write_cluster(fo, cluster_dir)

    def generate_cmd(self, mol:Molecule):
        return f"{self.executable} -b {self.inpfile} > {self.outfile}"

    def check_output(self, mol:Molecule):
        if not os.path.exists(self.result_frcmod):
            possibledir = os.path.basename(self.inpfile).replace(".inp", "")
            self.result_frcmod = os.path.join(self.result_dir, possibledir, "{}.frcmod".format(mol.name))
        if not os.path.exists(self.result_frcmod):
            print("Error: no result found for", mol.name)
            return False

    def process_output(self, mol:Molecule):
        mol.frcmod = self.result_frcmod
        mol.update()
        pass

    def run(self, mol:Molecule):
        self.check_system(mol)
        self.predict_output(mol)
        self.gen_forcefield(mol)
        self.generate_input(mol)
        cmd = self.generate_cmd(mol)
        print(cmd)
        self.execute(cmd)
        self.check_output(mol)
        self.process_output(mol)

