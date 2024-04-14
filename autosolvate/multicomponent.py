#--------------------------------------------------------------------------------------------------#
# multicomponent.py. 
# Description: 
#   This module can handle structure files containing multiple molecules.
#   1. Generates the lib and frcmod file for each separate molecule when it is not h2o or amino acids.
#   2. Generates the lib, prmtop and inpcrd files of this whole structure.
# Update 2022-02-04:
#   1. Can generate lib, mol2 and prmtop file for an xyz with mutiple fragments with charges
#   2. Can create solvent box for xyz file with multiple fragments.
#   3. Able to rearrange shuffled pdb files into ordered form that antechamber can process.
# author: Fangning Ren (2022-02-04) 
# path: autosolvate/multicomponent.py
#--------------------------------------------------------------------------------------------------#
import getopt, sys, os
import subprocess
from .molecule import *
from .dockers import *
from .utils import *

from autosolvate.autosolvate import *

class MulticomponentParamsBuilder():
    def __init__(self, 
                 xyzfile: str, name="", residue_name="SLU", charge=0, spinmult=1, 
                 charge_method="resp", folder = WORKING_DIR, **kwargs): 
        """
        Create amber parameter files for a single xyz or pdb file with multiple separate fragments

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any structural files that openbabel recognizes.
        name : array_like, Optional. default: the base name of the provided structure file.
            Not used
        residue_name : array_like, Optional. default: Residue name provided in pdb file or assigned as UAA, UAB, UAC, etc.
            Residue names for each fragments. A list of strings of three capital letters. Its length should equal to the number of fragments in xyzfile. If this parameter is not given, the residues will be assigned by "U" plus "AB","AC",..."AZ", "BA"...
        charge : dict | array_like, Optional. default: 0
            Charge for each fragment. A list of integer with its length equal to the number of fragments, or a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as neutral. 
        spinmult : dict | array_like, Optional. default: 0
            Multiplicity for each fragment. A list of integer with its length equal to the number of fragments, or a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as singlet. 
        outputFile : str, Optional, default='water_solvated'
            Filename-prefix for outputfiles
        pre_optimize_fragments : bool, Optional, default: False
            do geometry optimization with MMFF94 forcefield in OpebBabel before running antechamber
        srun_use : bool, Optional, default='False
            Run all commands with a srun prefix
        gaussianexe : str, Optional, default: g16
            name of the Gaussian executeble
        gaussiandir : str, Optional, default: $GAUSSIANDIR
            path of Gaussian
        amberhome : str, Optional, default: $AMBERHOME
            path of amber
        deletefiles : bool, Optional, default: False
            Delete all temporary files except the .prmtop and .inpcrd file of the pdb file provided.
        """
        
        self.folder = folder
        self.mol = MoleculeComplex(xyzfile, name=name, residue_name=residue_name, charges=charge, multiplicities=spinmult, folder = self.folder)
        self.charge_method = charge_method

        self.single_molecule_pipeline = [
            AntechamberDocker(charge_method = self.charge_method, workfolder = self.folder),
            ParmchkDocker(workfolder=self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.complex_pipeline = [TleapDocker(workfolder=self.folder)]

    def build(self):
        for m in self.mol.newmolecules:
            if self.charge_method == "resp":
                build_resp_terachem(m, folder = self.folder)
            else:
                for docker in self.single_molecule_pipeline:
                    docker.run(m)
        for docker in self.complex_pipeline:
            docker.run(self.mol)

       
class MulticomponentSolventBoxBuilder():
    def __init__(self, 
                 xyzfile:str, slu_charge=0, slu_spinmult=1, charge_method="resp", slu_count = 1,
                 solvent = "water", solvent_frcmod = "", solvent_off = "", solvent_box_name = "SLVBOX",
                 slv_generate = False, slv_xyz = "", slv_count = 210*8,
                 cube_size = 54, closeness = 0.8, folder = WORKING_DIR, outputFile = "",
                 **kwargs):

        self.kwargs = kwargs

        if not outputFile:
            outputFile = solvent + "_solvated"
        self.folder = folder
        if "slu_netcharge" in kwargs and isinstance(kwargs["slu_netcharge"], dict) and slu_charge == 0:
            slu_charge = kwargs["slu_netcharge"]
        self.solute = MoleculeComplex(xyzfile, slu_charge, slu_spinmult, folder = self.folder)
        self.charge_method = charge_method

        self.solvent = self.get_solvent(solvent, slv_xyz, solvent_frcmod, solvent_off, slv_generate, slv_count, solvent_box_name)
        self.system = SolvatedSystem(solvent + "_solvated", solute = self.solute, solvent=self.solvent,
                                     cubesize=cube_size, closeness=closeness, solute_number=slu_count, solvent_number=slv_count,
                                     folder = self.folder)

        self.single_molecule_pipeline = [
            AntechamberDocker(charge_method = self.charge_method, workfolder = self.folder),
            ParmchkDocker(workfolder=self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.complex_pipeline = [TleapDocker(workfolder=self.folder)]
        self.solvent_pipeline = [
            AntechamberDocker("bcc", "mol2", workfolder = self.folder),
            ParmchkDocker("frcmod", workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.prebuilt_solvation = [
            TleapDocker(workfolder = self.folder)
        ]
        self.custom_solvation = [
            PackmolDocker(workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]
    
    def get_solvent(self, solvent:str, slv_xyz:str = "", solvent_frcmod:str = "", solvent_off:str = "", slv_generate:bool = False, slv_count:int = 210*8, solvent_box_name:str = "SLVBOX"):
        if solvent in AMBER_SOLVENT_DICT:
            # amber solvents
            self_solvent = AMBER_SOLVENT_DICT[solvent]
        elif solvent in custom_solv_dict:
            # solvent data prepared by autosolvate
            solvPrefix = custom_solv_dict[solvent]
            solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".frcmod"))
            solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".prep"))
            solvent_pdb_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".pdb"))
            self_solvent = Molecule(solvent_pdb_path, 0, 1, solvent, residue_name = custom_solv_residue_name[solvent], folder = self.folder)
            self_solvent.frcmod = solvent_frcmod_path
            self_solvent.prep   = solvent_prep_path
            self_solvent.number = slv_count
        elif os.path.exists(solvent_frcmod) and os.path.exists(solvent_off):
            # using prebuilt solvent box
            self_solvent = SolventBox(solvent_off, solvent_frcmod, name = solvent, folder = self.folder, box_name=solvent_box_name)
        elif os.path.exists(slv_xyz) and slv_generate == True:
            # generate solvent box
            self_solvent = Molecule(slv_xyz, folder = self.folder)
            self_solvent.number = slv_count
        else:
            raise ValueError("Solvent not found")
        return self_solvent

    def build(self):
        for m in self.solute.newmolecules:
            if self.charge_method == "resp":
                build_resp_terachem(m, folder = self.folder)
            else:
                for docker in self.single_molecule_pipeline:
                    docker.run(m)
        for docker in self.complex_pipeline:
            docker.run(self.solute)

        if isinstance(self.solvent, Molecule):
            for docker in self.solvent_pipeline:
                docker.run(self.solvent)
        if isinstance(self.solvent, SolventBox):
            self.prebuilt_solvation[0].run(self.system)
        elif isinstance(self.solvent, Molecule):
            self.custom_solvation[0].run(self.system)
            self.custom_solvation[1].run(self.system)





if __name__ == "__main__":
    inst = MulticomponentParamsBuilder("PAHs.pdb", deletefiles=True)
    inst.build()
