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
from TleapDocker import TleapDocker 
import tools 


@dataclass 
class Molecule: 
    #constants 
    _SUPPORT_INPUT_FORMATS = ['pdb', 'xyz'] 
    
    #required arguments
    inputfile:    str
    name :        str  
    charge:       int
    multiplicity: int
    mol_type:     str

    #positional arguments 
    residue_name: str = 'MOL' 
    
    #set by __post_init__ 
    pdb:    str = field(init=False)
    
    #files 
    '''
    @TODO 
    make these thing store doc instead of a string
    or a specific place for storing files  
    '''
    mol2:   str = None
    frcmod: str = None
    lib:    str = None 
    prmtop: str = None
    inpcrd: str = None
    box:    str = None

    
    def __post_init__(self) -> None: 
        if self.inputfile == '': 
            raise Warning('inputfile is empty')

        if not isinstance(self.name, str): 
            raise Exception('name is not a string') 

        if self.charge not in [-1, 0, 1]: 
            raise Exception('invalid charge') 

        if self.multiplicity not in [1, 2, 3]: 
            raise Exception('invalid multiplicity')
            
        if self.mol_type not in ['solute', 'solvent', 'amber_solvent']:
            raise Exception('invalid mol_type')

        if self.mol_type != 'amber_solvent': 
            self.set_pdb()

    def set_output_folder(self) -> None:
        self.output_folder = WORKING_DIR + self.name + '/' 
        if not os.path.exists(self.output_folder):
            os.makedirs(self.output_folder)
        return 
        

    def set_pdb(self) -> None:
        basename, name, ext = tools.extract_basename_name_extension(self.inputfile)
        if ext not in self._SUPPORT_INPUT_FORMATS: 
            raise Exception('Input file format not supported in class Molecule()') 
        if ext == 'pdb': 
            self.pdb = self.inputfile 
        if ext == 'xyz': 
            self.pdb = tools.convert_xyz_to_pdb(self.inputfile) 

    
    def update(self) -> None:
        workspace = WORKING_DIR + '{}/**/*'.format(self.name) 
        
        if self.name + '.mol2' in glob.glob(
            workspace + '.mol2', recursive=True
        ):
            self.mol2 = self.name + '.mol2' 

        if self.name + '.frcmod' in glob.glob(
            workspace + '.frcmod', recursive=True
        ):
            self.frcmod = self.name + '.frcmod' 
        
        if self.residue_name + '.lib' in glob.glob(
            workspace + '.lib', recursive=True
        ):
            self.lib = self.residue_name + '.lib' 

        # if self.name + '.prmtop' in glob.glob(workspace + '.prmtop', recursive=True): 
        #     self.prmtop = self.name + '.prmtop' 

        # if self.name + '.inpcrd' in glob.glob(workspace + '.inpcrd', recursive=True): 
        #     self.inpcrd = self.name + '.inpcrd' 

        

def update_mol(mol: object) -> None:
    r'''
    @TODO: 
        1. implement update for solute, or there is no need to update solute 
    '''
    if mol.mol_type == 'solvent': 
        update_solvent(mol)
    else: 
        print('not implemented yet')
    return 
        


def update_solvent(mol: object) -> None:
    if mol.mol2 == None or mol.frcmod == None:
        AntechamberDocker().run(mol)
    if mol.lib == None: 
        TleapDocker().run(mol)
    mol.update()
    return 

