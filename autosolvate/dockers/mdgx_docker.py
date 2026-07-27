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

class MDGXDocker(GeneralDocker):
    def __init__(self, workfolder:str="wlp"):
        super().__init__(executable="mdgx", workfolder=workfolder)
        
        os.makedirs(self.workfolder, exist_ok=True)

        self.target_frcmod  = ""
        self.target_mol2    = ""
        self.result_frcmod  = ""

        self.energys = []
        self.conformers = []

        self.qmdatafile = os.path.join(self.workfolder, "qmdata.dat")
        self.mmtrajfile = os.path.join(self.workfolder, "mmtraj.netcdf")
        self.accrepfile = os.path.join(self.workfolder, "fit.m")
        self.inpfile = os.path.join(self.workfolder, "inpfile.inp")
        self.outfile = os.path.join(self.workfolder, "outfile.out")

    def add_data(self, mol:Molecule, qmfolder:str, prefix:str = ""):
        prefix = prefix if prefix else mol.name
        xyzs, grds, enes = read_qm_folder(qmfolder, prefix)
        nconf = len(xyzs)
        if nconf == 0:
            print("No data found for", mol.name)
            return
        self.energys += enes
        self.conformers += xyzs
    
    def write_qdata(self, mol:Molecule):
        with open(self.qmdatafile, "w") as f:
            f.write( "% The following energies apply to conformations described by the topology \n")
            f.write(f"% {mol.prmtop}.  Structures are given in trajectory {self.mmtrajfile}, count: {len(self.energys)}.\n")
            f.write( "\n")
            f.write( "% Item 0:\n")
            f.write(f"%   Energy -> {mol.name}*.out\n")
            for e in self.energys:
                f.write("{:>16.8f}\n".format(e))

    def write_netcdf(self, mol:Molecule):
        t0 = md.load_xyz(mol.xyz, top = mol.prmtop)
        cfms = []
        for xyz in self.conformers:
            t = deepcopy(t0)
            t.xyz[0] = xyz / 10.0
            cfms.append(t)
        traj = md.join(cfms)
        traj.save_netcdf(self.mmtrajfile)

    def check_system(self, mol:Molecule):
        pass

    def copy_frcmod(self, mol:Molecule):
        # copy the molecular frcmod file to the workfolder
        shutil.copy(mol.frcmod, self.workfolder)
        os.rename(os.path.join(self.workfolder, os.path.basename(mol.frcmod)), self.target_frcmod)

    def predict_output(self, mol:Molecule):
        self.inpfile            = os.path.join(self.workfolder, "{}-mdgxfit.inp"    .format(mol.name))
        self.outfile            = os.path.join(self.workfolder, "{}-mdgxfit.out"    .format(mol.name))
        self.qmdatafile         = os.path.join(self.workfolder, "{}-qmdata.dat"     .format(mol.name))
        self.mmtrajfile         = os.path.join(self.workfolder, "{}-mmtraj.netcdf"  .format(mol.name))
        self.accrepfile         = os.path.join(self.workfolder, "{}-fit.m"          .format(mol.name))
        self.target_frcmod      = os.path.join(self.workfolder, "{}-old.frcmod"     .format(mol.name))
        self.result_frcmod      = os.path.join(self.workfolder, "{}-new.frcmod"     .format(mol.name))

    def write_input_head(self, mol:Molecule, fo:TextIO):
        fo.write("&files\n")
        fo.write("  -parm /home/fren5/psi4conda/envs/autosolvate/dat/leap/parm/gaff.dat\n")
        fo.write("  -fmod {}\n".format(self.target_frcmod))
        fo.write("  -d    {}\n".format(self.result_frcmod))
        fo.write("  -o    {}\n".format(self.outfile))
        fo.write("&end\n")
        fo.write("\n")

    def write_input_bond_fit(self, mol:Molecule, fo:TextIO, bkeys:List[str]):
        fo.write("  % Bond fit\n")
        for bkey in bkeys:
            fo.write("  fitb        {:<2s} {:<2s}\n".format(bkey[0], bkey[1]))
        fo.write("  FitBondEq   1,       \n")
        fo.write("  brst        0.0002,  \n")
        fo.write("  brstcpl     1.0,     \n")

    def write_input_angle_fit(self, mol:Molecule, fo:TextIO, akeys:List[str]):
        fo.write("  % Angle fit\n")
        for akey in akeys:
            fo.write("  fita        {:<2s} {:<2s} {:<2s}\n".format(akey[0], akey[1], akey[2]))
        fo.write("  FitAnglEq   1,      \n")
        fo.write("  arst        0.0002,  \n")
        fo.write("  arstcpl     1.0,     \n")

    def write_input_torsion_fit(self, mol:Molecule, fo:TextIO, tkeys:List[str]):
        fo.write("  % Torsion fit\n")
        for tkey in tkeys:
            fo.write("  fith        {:<2s} {:<2s} {:<2s} {:<2s}\n".format(tkey[0], tkey[1], tkey[2], tkey[3]))
        fo.write("  fitD                \n")
        fo.write("  hrst        0.0002,  \n")
    
    def write_input_general(self, mol:Molecule, fo:TextIO):
        fo.write("&param\n")
        fo.write("  System {} {} {}     \n".format(mol.prmtop, self.mmtrajfile, self.qmdatafile))
        fo.write("  ParmOutput frcmod   \n")
        fo.write("  eunits hartree,     \n")
        fo.write("  accrep {}           \n".format(self.accrepfile))
        fo.write("  verbose 2,          \n")
        fo.write("\n")

    def generate_input(self, mol:Molecule):
        self.write_qdata(mol)
        self.write_netcdf(mol)

        fmod = FrcmodFile(self.target_frcmod)
        bkeys, akeys, dkeys, ikeys = fmod.get_mdgx_labels()
        with open(self.inpfile, "w") as fo:
            self.write_input_head(mol, fo)
            self.write_input_general(mol, fo)
            self.write_input_bond_fit(mol, fo, bkeys)
            self.write_input_angle_fit(mol, fo, akeys)
            self.write_input_torsion_fit(mol, fo, dkeys)
            fo.write("&end\n")

    def generate_cmd(self, mol:Molecule):
        return f"{self.executable} -i {self.inpfile}"

    def check_output(self, mol:Molecule):
        pass

    def process_output(self, mol:Molecule):
        oldfmod = FrcmodFile(self.target_frcmod)
        newfmod = FrcmodFile(self.result_frcmod)
        oldfmod.update(newfmod)
        oldfmod.write(mol.frcmod)
        mol.update()

    def run(self, mol:Molecule):
        self.check_system(mol)
        self.predict_output(mol)
        self.copy_frcmod(mol)
        self.generate_input(mol)
        cmd = self.generate_cmd(mol)
        print(cmd)
        self.execute(cmd)
        self.check_output(mol)
        self.process_output(mol)
