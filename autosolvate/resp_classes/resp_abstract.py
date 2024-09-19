"""
resp.py

Generate resp charge with various quantum chemistry packages

Handles the primary functions
"""
from abc import ABC, abstractmethod
from openbabel import openbabel as ob
from openbabel import pybel
import os
import re
import copy
from autosolvate.globs import keywords_avail, available_qm_programs


class RespABC(ABC):
    ### Attributes that are in common for different QM packages should go in here.
    def __init__(self, **kwargs):
        print("inputs dictionary:", kwargs)
        self.keywords_avail = keywords_avail
        self.pdbfile = False
        self.molname = kwargs["molname"] if "molname" in kwargs.keys() else "undef"
        self.qm_program = kwargs["qm_program"] if "qm_program" in kwargs.keys() else "gamess"
        self.qm_exe = kwargs["qm_exe"] if "qm_exe" in kwargs.keys() else None
        self.qm_dir = kwargs["qm_dir"] if "qm_dir" in kwargs.keys() else None
        self.rundir = kwargs["rundir"] if "rundir" in kwargs.keys() else os.getcwd()
        self.rundir = os.path.abspath(self.rundir)
        self.resp_scr_dir = os.path.join(os.path.abspath(self.rundir), 'resp_scr')
        self.initialization_check(**kwargs)
        self.xyzfile = kwargs["xyzfile"] if "xyzfile" in kwargs.keys() else "undef"
        self.nprocs = kwargs["nprocs"] if "nprocs" in kwargs.keys() else "undef"

    def initialization_check(self, **kwargs):
        for key, value in kwargs.items():
            if key not in self.keywords_avail:
                raise KeyError("RESP charge fitting: unrecoganized key")

        if self.qm_program not in available_qm_programs:
            raise ValueError("RESP charge fitting: specified program is not supported yet")

        if "pdbfile" in kwargs.keys():
            self.pdbfile = kwargs["pdbfile"]
            self.charge = kwargs["charge"] if "charge" in kwargs.keys() else 0
            self.spinmult = kwargs["spinmult"] if "spinmult" in kwargs.keys() else 1
            self.molecule = self.pdb2obmol(self.pdbfile, charge=self.charge, spinmult=self.spinmult)
            if self.molname == "undef":
                  self.molname = self.pdb2molname(self.pdbfile) 
        else:
            raise KeyError("Error: Must input a pdbfile file for resp charge fitting.")

    def pdb2obmol(self, pdbfile, charge=0, spinmult=1):
        molecule = pybel.readfile('pdb', pdbfile).__next__()
        obmol = molecule.OBMol
        if charge != 0:
           obmol.SetTotalCharge(charge)
        if spinmult > 1:
           obmol.SetTotalSpinMultiplicity(spinmult)
        return obmol
     
    def pdb2molname(self, pdbfile):
        if pdbfile.find('.pdb') == -1:
             raise ValueError('The input argument must be a file name ended with .pdb')
        molname = copy.copy(pdbfile)
        molname = molname.replace('.pdb','')
        molname = molname.replace('.xyz','')
        return molname

