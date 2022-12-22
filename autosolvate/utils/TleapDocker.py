from openbabel import pybel
from openbabel import openbabel as ob
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools
from Molecule import Molecule 

class TleapDocker: 
    _ff14SB  = 'leaprc.ff14SB'               #Source leaprc file for ff14SB protein force field
    _gaff    = 'leaprc.gaff'                 #Source leaprc file for gaff force field


    def __init__(self):
        pass                        


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
        cmd = 'tleap -s -f leap_{}.in > leap_{}.out'.format(mol.name, mol.name)
        return cmd


    def write_tleap_in(self, mol: Molecule) -> None:  
        f = open('leap_{}.in'.format(mol.name), 'w')
        self.add_forcefield(f)
        self.add_frcmod(f, mol)
        self.add_mol2(f, mol)
        self.write_set_head_and_tail(f)
        self.write_save_outputs(f, mol)
        self.write_end(f)
        f.close() 
        pass 


    def add_forcefield(self, doc: object) -> None:
        doc.write('{:<20}  {:<20} \n'.format('source', self._ff14SB)) 
        doc.write('{:<20}  {:<20} \n'.format('source', self._gaff)) 


    def add_frcmod(self, doc: object, mol: Molecule) -> None: 
        doc.write('{:<20}  {:<20} \n'.format('loadamberparams', mol.frcmod))


    def add_mol2(self, doc: object, mol: Molecule) -> None:
        doc.write('{:<5} {:<15} {:<20} \n'.format(mol.residue_name, '= loadmol2', mol.mol2))
    

    def write_set_head_and_tail(self, doc: object) -> None:
        print('set head and tail is not implemented yet') 
        pass 


    def write_save_outputs(self, doc: object, mol: Molecule) -> None: 
        doc.write('{:<20} {:<5} {:<15} \n'.format(
            'saveoff', mol.residue_name, mol.name + '.lib')
        )
        doc.write('{:<20} {:<5} {:<15} \n'.format(
            'savepdb', mol.residue_name, mol.name + '.pdb')
        ) 
        doc.write('{:<20} {:<5} {:<15} {:<15} \n'.format(
            'saveamberparm', mol.residue_name, mol.name + '.prmtop', mol.name + '.inpcrd')
        )


    def write_end(self, doc):
        doc.write('quit') 


    def check_mol(self, mol: Molecule) -> None: 
        if mol.name is None: 
            raise Exception('name is not set') 
        if mol.mol2 is None:
            raise Exception('mol2 file is not set')
        if mol.frcmod is None:
            raise Exception('frcmod file is not set')
        if mol.residue_name is None:
            raise Exception('residue name is not set')
        pass