from openbabel import pybel
from openbabel import openbabel as ob
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools


class ParmchkDocker:
    def __init__(self, 
                 out_format: str = 'frcmod'
    ) -> None:
        self.out_format = out_format


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
        $AMBERHOME/bin/parmchk2 -i 1.prmtop -f mol2 -o 1.frcmod
        '''
        cmd =  self.set_executable()    + ' ' 
        cmd += self.set_input(mol)      + ' ' 
        cmd += self.set_output(mol)     + ' ' 
        return cmd 


    def set_executable(self) -> str:
        return '$AMBERHOME/bin/parmchk2'


    def set_input(self, mol: object) -> str:    
        if mol.mol2 is not None:
            return '-i %s -f %s' % (os.path.basename(mol.mol2), 'mol2') 

        raise Exception('only support mol2 as input format') 


    def set_output(self, mol: object) -> str: 
        return '-o %s.%s' % (mol.name, self.out_format)

       


def check_mol(*args) -> None:
    for mol in args:
        if mol.name is None:
            raise Exception('mol.name is None')