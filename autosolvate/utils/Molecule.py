# @TODO
# 1. modify the @EXAMPLE in the docstring 
# 2. implement charge and multiplicity determination
# 3. TleapDocker().run(mol) not fully implemented yet 
# 4. check update() method in Molecule class 

from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np
import getopt, sys, os, subprocess, pkg_resources, glob 
from dataclasses import dataclass, field, asdict
from Common import * 
from AntechamberDocker import AntechamberDocker 
import tools 


@dataclass 
class Molecule: 
    #constants 
    SUPPORT_INPUT_FORMATS = ['pdb', 'xyz'] 
    
    #required arguments
    inputfile:    str 
    charge:       int 
    multiplicity: int 
    
    #positional arguments 
    residue_name: str = 'MOL' 
    Is_a_solvent: bool = False 
    
    #set by __post_init__ 
    pdb:    str = field(init=False)
    name:   str = field(init=False) 
    
    #files 
    mol2:   str = None
    frcmod: str = None
    lib:    str = None 
    prmtop: str = None
    inpcrd: str = None

    
    def __post_init__(self):
        self.set_pdb() 
        self.set_name() 


    def set_pdb(self) -> None:
        basename, name, ext = tools.extract_basename_name_extension(self.inputfile)
        if ext not in self.SUPPORT_INPUT_FORMATS: 
            raise Exception('Input file format not supported in class Molecule()') 
        if ext == 'pdb': 
            self.pdb = self.inputfile 
        if ext == 'xyz': 
            self.pdb = tools.convert_xyz_to_pdb(self.inputfile) 


    def set_name(self) -> None: 
        basename, name, ext = tools.extract_basename_name_extension(self.inputfile)
        if name == '': 
            raise Exception('invalid input file name') 
        self.name = name 

    
    def update(self) -> None:
        current_dir = os.getcwd() + '/**/*'
        
        if self.name + '.mol2' in glob.glob(
            current_dir + '.mol2', recursive=True
        ):
            self.mol2 = self.name + '.mol2' 

        if self.name + '.frcmod' in glob.glob(
            current_dir + '.frcmod', recursive=True
        ):
            self.frcmod = self.name + '.frcmod' 
        
        if self.residue_name + '.lib' in glob.glob(
            current_dir + '.lib', recursive=True
        ):
            self.lib = self.residue_name + '.lib' 

        # if self.name + '.prmtop' in glob.glob(current_dir + '.prmtop', recursive=True): 
        #     self.prmtop = self.name + '.prmtop' 

        # if self.name + '.inpcrd' in glob.glob(current_dir + '.inpcrd', recursive=True): 
        #     self.inpcrd = self.name + '.inpcrd' 

        


def update_mol(mol: Molecule) -> Molecule:
    if mol.Is_a_solvent == True: 
        return update_solvent(mol) 


def update_solvent(mol: Molecule) -> Molecule:
    if mol.mol2 == None or mol.frcmod == None:
        AntechamberDocker().run(mol)
    if mol.lib == None: 
        TleapDocker().run(mol)
    return mol


