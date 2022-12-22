# @TODO: 
# 1. change input to a Molecule object, remove charge and multiplicity 
# 2. this does not need to be a class 
# 3. check if 'mol: Molecule' usage is valid 

from openbabel import pybel
from openbabel import openbabel as ob
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools
from Molecule import Molecule  

class AntechamberDocker:

    SUPPORT_INPUT_FORMATS           = ['pdb']
    SUPPORT_CHARGE_FITTING_METHODS  = ['bcc']

    def __init__(self, charge_fiiting_method:str = 'bcc', output_format:str = 'mol2') -> None:      
        #setting
        self.output_format                  = output_format
        self.charge_fiiting_method          = charge_fiiting_method  


    def run(self, mol: Molecule):
        self.check_mol(mol) 
        cmd = self.generate_cmd(mol)
        if DRY_RUN:
            print(cmd) 
            return
        else:
            subprocess.call(cmd, shell=True)     


    @tools.srun()
    def generate_cmd(self, mol: Molecule) -> str:
        '''
        @EXAMPLE: 
        $AMBERHOME/bin/antechamber -i 1.pdb -fi pdb -o 1.mol2 -fo mol2 -c bcc -nc 0 -m 1 -rn MOL
        '''
        cmd =  self.set_executable(mol)    + ' '
        cmd += self.set_input(mol)         + ' '
        cmd += self.set_output(mol)        + ' ' 
        cmd += self.set_charge(mol)        + ' ' 
        cmd += self.set_multiplicity(mol)  + ' ' 
        cmd += self.set_charge_method(mol) + ' ' 
        cmd += self.set_residue_name(mol)  + ' ' 
        return cmd 


    def set_executable(self) -> str :
        return '$AMBERHOME/bin/antechamber'
        

    def set_input(self, mol: Molecule) -> str: 
        pdb = mol.pdb 
        return '-i %s -fi %s' % (pdb, 'pdb') 


    def set_output(self, mol: Molecule) -> str:
        return '-o %s.mol2 -fo mol2' % mol.name 


    def set_charge(self, mol: Molecule) -> str:
        if mol.charge == 0:
            return '' 
        return '-nc %d' % mol.charge 


    def set_multiplicity(self, mol: Molecule) -> str: 
        if mol.multiplicity == 1:
            return '' 
        return '-m %d' % mol.multiplicity 
        

    def set_charge_method(self) -> str: 
        if self.charge_method not in self.SUPPORT_CHARGE_FITTING_METHODS:
            raise Exception('Charge fitting method not supported in class AntechamberDocker()')
        return '-c %s' % self.charge_method


    def set_residue_name(self, mol: Molecule) -> str: 
        if self.residue_name == 'MOL':
            return '' 
        return '-rn %s' % mol.residue_name


    def check_mol(self, mol: Molecule) -> None:
        if mol.pdb is None:
            raise Exception('mol.pdb is None')
        if mol.name is None: 
            raise Exception('mol.name is None')
        if mol.charge is None: 
            raise Exception('mol.charge is None')
        if mol.multiplicity is None: 
            raise Exception('mol.multiplicity is None')
        if mol.residue_name is None: 
            raise Exception('mol.residue_name is None')
        


if __name__ == '__main__': 
    import doctest

    global DRY_RUN 
    DRY_RUN = True

    doctest.testmod() 
    