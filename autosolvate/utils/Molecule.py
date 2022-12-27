# @TODO
# 1. modify the @EXAMPLE in the docstring 
# 2. implement charge and multiplicity determination
# 3. TleapDocker().run(mol) not fully implemented yet 
# 4. check update() method in Molecule class 
# @DISCUSSION 
# 1. should we store doc object instead of a string of the filename?
#    for mol2, frcmod, lib, prmtop, inpcrd, box
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
    name :          str  
    charge:         int
    multiplicity:   int
    mol_type:       str

    #positional arguments 
    residue_name:   str = 'MOL'
    #files 
    pdb:            str = None
    xyz:            str = None
    mol2:           str = None
    frcmod:         str = None
    lib:            str = None 
    prmtop:         str = None
    inpcrd:         str = None
    box:            str = None

   
    def __post_init__(self) -> None:
        ''' 
        @NOTE:
        1. amber dic solvent can not use os.path.exists() to check if the frcmod exists 
        
        2. some frcmod files is named as frcmod.*  (not *.frcmod)

        3. remember molecule can be initailzed many ways: 
            1. from pdb file
            2. from xyz file (not implemented yet) 
            3. from amber library. (we provide lib, frcmod and box files) 
            4. from mol2 and frcmod files (we provide these files in data/) 

        @TODO:
        1. custome solvent does have mol2 and frcmod files.  
           we need to check if mol2 and frcmod files provided in the working directory or data/ 
        '''

        if not isinstance(self.name, str): 
            raise Exception('name is not a string') 

        if self.charge not in [-1, 0, 1]: 
            raise Exception('invalid charge') 

        if self.multiplicity not in [1, 2, 3]: 
            raise Exception('invalid multiplicity')
            
        if self.mol_type not in ['solute', 'solvent']:
            raise Exception('invalid mol_type')

        if not isinstance(self.residue_name, str): 
            raise Exception('residue_name is not a string') 

        if self.pdb is not None:
            if not isinstance(self.pdb, str): 
                raise Exception('pdb is not a string') 
            if not os.path.exists(self.pdb): 
                raise Exception('pdb file does not exist') 

        if self.xyz is not None: 
            if not isinstance(self.xyz, str): 
                raise Exception('xyz is not a string') 
            if not os.path.exists(self.xyz): 
                raise Exception('xyz file does not exist') 

        if self.mol2 is not None: 
            if not isinstance(self.mol2, str): 
                raise Exception('mol2 is not a string') 
                
        if self.frcmod is not None:
            if not isinstance(self.frcmod, str): 
                raise Exception('frcmod is not a string')

        if self.lib is not None:
            if not isinstance(self.lib, str): 
                raise Exception('lib is not a string') 

        if self.prmtop is not None:
            if not isinstance(self.prmtop, str): 
                raise Exception('prmtop is not a string') 
            
        if self.inpcrd is not None: 
            if not isinstance(self.inpcrd, str): 
                raise Exception('inpcrd is not a string') 

        if self.box is not None: 
            if not isinstance(self.box, str): 
                raise Exception('box is not a string') 


    def __eq__(self, o: object) -> bool: 
        ''' 
        @TODO 
        test this method 
        ''' 
        if not isinstance(o, Molecule):
            return False 
        for key in asdict(self).keys(): 
            if getattr(self, key) != getattr(o, key): 
                return False
        return True


    # def set_output_folder(self) -> None:
    #     self.output_folder = WORKING_DIR + self.name + '/' 
    #     if not os.path.exists(self.output_folder):
    #         os.makedirs(self.output_folder)
    #     return 
        
    def update(self) -> None:
        search_range = WORKING_DIR + '{}/**/*'.format(self.name) 
        
        if self.name + '.mol2' in glob.glob(
            search_range + '.mol2', recursive=True
        ):
            self.mol2 = self.name + '.mol2' 

        if self.name + '.frcmod' in glob.glob(
            search_range + '.frcmod', recursive=True
        ):
            self.frcmod = self.name + '.frcmod' 
        
        if self.residue_name + '.lib' in glob.glob(
            search_range + '.lib', recursive=True
        ):
            self.lib = self.residue_name + '.lib' 

        # if self.name + '.prmtop' in glob.glob(search_range + '.prmtop', recursive=True): 
        #     self.prmtop = self.name + '.prmtop' 

        # if self.name + '.inpcrd' in glob.glob(search_range + '.inpcrd', recursive=True): 
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

