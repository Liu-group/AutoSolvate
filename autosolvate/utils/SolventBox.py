# @TODO
# 1. modify the @EXAMPLE in the docstring 
# 2. implement charge and multiplicity determination
# 4. check update() method in Molecule class 

# @DISCUSSION 
# 1. should we store doc object instead of a string of the filename?
#    for mol2, frcmod, lib, prmtop, inpcrd, box

from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np
import getopt, sys, os, subprocess, pkg_resources, glob, shutil
from dataclasses import dataclass, field, asdict
from datetime import datetime
from Common import * 
from Molecule import Molecule, AMBER_SOLVENT_LIST 
import tools 

@dataclass 
class SolventBox:
    solute_list:            list  = field(default_factory=list, init=False)
    solvent_list:           list  = field(default_factory=list, init=False)

    cubesize:               int   = 54 
    closeness:              float = 0.8  

    duplicate_solute_num:   int   = 1  
    duplicate_solvent_num:  int   = 1 

    system_pdb:             str   = field(default=None, init=False)   
    name:                   str   = None
    folder:                 str   = field(default=None, init=False) 


    def __post_init__(self) -> None:
        self.set_box_name()
        self.set_folder()

        
    def set_box_name(self) -> str:
        '''
        @TODO: 
        name can be improved  
        but I have no idea now what naming scheme should be 
        '''
        if self.name is None: 
            now = datetime.now()
            current_time = now.strftime("%H_%M_%S")
            self.name = 'box_' + current_time


    def set_folder(self) -> None:
        path = WORKING_DIR + self.name + '/'
        if not os.path.exists(path):
            os.mkdir(path)
        self.folder = path 


    def add_solute(self, mol: object) -> None: 
        '''
        @TODO 
        check solute validity  
        ''' 
        self.solute_list.append(mol)
        self.update()


    def add_solvent(self, mol: object) -> None: 
        '''
        @TODO 
        check solvent validity 
        ''' 
        self.solvent_list.append(mol)
        self.update()


    def set_closeness(self, solvent: object, automate: bool = False) -> None: 
        r'''
        @TODO: 
            1. It only supports one solvent 
        '''     
        if automate: 
            if solvent.name == 'acetonitrile':
                self.closeness = 1.88 
            elif solvent.name == 'water': 
                self.closeness = 0.50
            elif solvent.name == 'methanol': 
                self.closeness = 0.60 
            elif solvent.name == 'nma': 
                self.closeness = 0.58 
            elif solvent.name == 'chloroform': 
                self.closeness = 0.58 
            else: 
                raise Warning('unknown solvent name')
        else: 
            pass 

    
    def check_residue_name(self) -> None: 
        '''
        @TODO 
        check residue name 
        ''' 
        pass 


    def update(self) -> None:
        '''
        @TODO: 
        simplify this code
        '''
        search_range = WORKING_DIR + '/**/*'
        
        for file in glob.glob(search_range, recursive=True):
            
            if file.endswith('system.pdb'):
                self.system_pdb = system.pdb
                if self.name+'/' not in file:
                    shutil.copy(file, WORKING_DIR + self.name + '/')
            
            for solute in self.solute_list:
                if file.endswith(solute.name + '.pdb'):
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')
                if file.endswith(solute.name + '.mol2'):
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                if file.endswith(solute.name + '.frcmod'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                if file.endswith(solute.name + '.lib'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')
                if file.endswith(solute.name + '.prmtop'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                if file.endswith(solute.name + '.inpcrd'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 

            for solvent in self.solvent_list: 
                if file.endswith(solvent.name + '.pdb'):
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')
                if file.endswith(solvent.name + '.mol2'):
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                if file.endswith(solvent.name + '.frcmod'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                if file.endswith(solvent.name + '.lib'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')
                if file.endswith(solvent.name + '.prmtop'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                if file.endswith(solvent.name + '.inpcrd'): 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')