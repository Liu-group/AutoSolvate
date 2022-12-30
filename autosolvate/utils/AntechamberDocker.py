# @TODO: 
# 1. do we want output format to be .prep? 
# 2. add parmchk command (maybe a class)

# @NOTE: 
# 1. check if 'mol: object' usage is valid 
#    - yes, I think it is valid 
# 2. this does not need to be a class 
#    - yes it is a functional oriented class 
from openbabel import pybel
from openbabel import openbabel as ob
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools


class AntechamberDocker:

    _SUPPORT_INPUT_FORMATS           = ['pdb']
    _SUPPORT_CHARGE_FITTING_METHODS  = ['bcc']

    def __init__(self, 
                 charge_fiiting_method: str = 'bcc', 
                 out_format:            str = 'mol2'
    ) -> None:      
        #setting
        self.out_format                 = out_format
        self.charge_fiiting_method      = charge_fiiting_method  


    def run(self, mol: object) -> None:
        check_mol(mol)    
        os.chdir(mol.name)
        cmd = self.generate_cmd(mol)
        if DRY_RUN:
            print(cmd) 
            return
        else:
            subprocess.call(cmd, shell=True)    
        os.chdir(WORKING_DIR)  


    @tools.srun()
    def generate_cmd(self, mol: object) -> str:
        '''
        @EXAMPLE: 
        $AMBERHOME/bin/antechamber -i 1.pdb -fi pdb -o 1.mol2 -fo mol2 -c bcc -nc 0 -m 1 -rn MOL
        '''
        cmd =  self.set_executable()       + ' '
        cmd += self.set_input(mol)         + ' '
        cmd += self.set_output(mol)        + ' ' 
        cmd += self.set_charge(mol)        + ' ' 
        cmd += self.set_multiplicity(mol)  + ' ' 
        cmd += self.set_charge_method() + ' ' 
        cmd += self.set_residue_name(mol)  + ' ' 
        print(cmd) 
        return cmd 


    def set_executable(self) -> str :
        return '$AMBERHOME/bin/antechamber'
        

    def set_input(self, mol: object) -> str:
        return '-i %s -fi %s' % (os.path.basename(mol.pdb), 'pdb') 

 
    def set_output(self, mol: object) -> str:
        '''
        @NOTE: 
        dont use mol.mol2 because it is None at this point 
        '''
        return '-o %s.%s -fo %s' % (mol.name, self.out_format, self.out_format) 
       

    def set_charge(self, mol: object) -> str:
        if mol.charge == 0:
            return '' 
        return '-nc %d' % mol.charge 


    def set_multiplicity(self, mol: object) -> str: 
        if mol.multiplicity == 1:
            return '' 
        return '-m %d' % mol.multiplicity 
        

    def set_charge_method(self) -> str: 
        if self.charge_fiiting_method not in self._SUPPORT_CHARGE_FITTING_METHODS:
            raise Exception('Charge fitting method not supported in class AntechamberDocker()')
        return '-c %s' % self.charge_fiiting_method


    def set_residue_name(self, mol: object) -> str: 
        if mol.residue_name == 'MOL':
            return '' 
        return '-rn %s' % mol.residue_name



def check_mol(mol: object) -> None:
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
    if os.path.exists(mol.name) is False:
        raise Exception('directory %s does not exist' % mol.name) 
        



if __name__ == '__main__': 
    import doctest

    global DRY_RUN 
    DRY_RUN = True

    doctest.testmod() 
    