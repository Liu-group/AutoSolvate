from openbabel import pybel
from openbabel import openbabel as ob
from multipledispatch import dispatch 
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools
from Molecule import Molecule 

class TleapDocker: 

    _ff14SB         = 'leaprc.ff14SB'               #Source leaprc file for ff14SB protein force field
    _gaff           = 'leaprc.gaff'                 #Source leaprc file for gaff force field
    _water_tip3p    = 'leaprc.water.tip3p'


    def __init__(self):     
        '''
        @NOTE:
        1. do not store any molecule obejct in this class attribute 
           because this class is a function class, not a data class 
        '''
        self._solute_num    = 0
        self._solvent_num   = 0 
        self._load_ions     = False
        

    def submit(self, cmd: str) -> None: 
        if DRY_RUN:
            print(cmd) 
            return
        else:
            subprocess.call(cmd, shell=True)


    @dispatch(object)
    def run(self, mol: object) -> None:
        check_args(mol)
        if self._solute_num == 1:
            self.write_tleap_in(mol)
            self.generate_cmd(mol)
        else: 
            raise Exception('invalid number of solute')
        self.submit(cmd)        

    @dispatch(object, object) 
    def run(self, mol1: object, mol2: object) -> None:
        r'''
        @TODO 
        1. finish this function 
        2. need to add to check if mol1 and mol2 are solute and solvent 
        3. need to check if solvent is custom or not 
        '''
        check_args(mol1, mol2)


    def check_args(self, *args: object) -> None: 
        for mol in args:
            check_mol_arributes(mol)                
        self._solute_num    = count_solute(*args)
        self._solvent_num   = count_solvent(*args)


    @tools.srun()
    def generate_cmd(self, mol: object) -> str:
        cmd = 'tleap -s -f leap_{}.in > leap_{}.out'.format(mol.name, mol.name)
        return cmd


    @dispatch(object)
    def write_tleap_in(self, mol: object) -> None:  
        r'''
        @QUESTION 
        1. Is check mol necessary?, what is it for? 
        '''
        f = open('leap_{}.in'.format(mol.name), 'w')
        self.add_forcefield(f)
        self.load_mol(f, mol, frcmod=True, mol2=True)  
        f.write('{:<20} {:<5} {:<15}        \n'.format('saveoff', mol.residue_name, mol.name + '.lib'))
        f.write('{:<20} {:<5} {:<15}        \n'.format('savepdb', mol.residue_name, mol.name + '.pdb'))
        f.write('{:<20} {:<5} {:<15} {:<15} \n'.format('saveamberparm', mol.residue_name, mol.name + '.prmtop', mol.name + '.inpcrd')) 
        f.write('{:<20}                     \n'.format('quit'))
        f.close() 

    @dispatch(object, object)
    def write_tleap_in(self, solute: object, solvent: object) -> None:
        r'''
        @param slu: solute MOLECULE object
        @param slv: solvent MOLECULE object

        @QUESTION: 
        1. does solvent need mol2 file? 
        '''
        #setting
        self._load_ions = True 
        #write tleap.in 
        f = open('leap_{}.in'.format('solventbox'), 'w')
        self.add_forcefield(f) 
        self.load_mol(f, solvent, frcmod=True, lib=True)
        self.load_mol(f, solute,  frcmod=True, mol2=True, check_mol=True)
        

    def add_forcefield(self, doc: object) -> None:
        doc.write('{:<20}  {:<20}       \n'.format('source', self._ff14SB)) 
        doc.write('{:<20}  {:<20}       \n'.format('source', self._gaff))
        if self._load_ions:
            doc.write('{:<20}  {:<20}   \n'.format('source', self._water_tip3p))
        

    def load_mol(self,  doc: object, mol: object, 
                        frcmod: bool = False, mol2: bool = False, lib: bool = False, check_mol: bool = False
    ) -> None: 
        r'''
        @QUESTION:
        1. what is check mol for? 
        '''
        if frcmod: 
            doc.write('{:<20}  {:<20}       \n'.format('loadamberparams', mol.frcmod))       
        if mol2:  
            doc.write('{:<5} {:<15} {:<20}  \n'.format(mol.residue_name, '= loadmol2', mol.mol2))
        if lib: 
            doc.write('{:<20}  {:<20}       \n'.format('loadoff', mol.lib)) 
        if check_mol: 
            doc.write('{:<20}  {:<20}       \n'.format('check', 'mol'))

    
    


def check_mol_arributes(mol: object) -> None: 
    '''
    @Description: 
        check if important attributes of the molecule are set for tleap 
    ''' 
    if mol.residue_name is None:
        raise Exception('residue name is not set')     
    if mol.mol2 is None:
        raise Exception('mol2 file is not set')
    if mol.frcmod is None:
        raise Exception('frcmod file is not set')


def count_solvent(*args: object) -> int: 
    ''' 
    @Description: 
        return the number of solvent molecules 
    '''
    solvent_num = 0 
    for mol in args:
        if mol.mol_type == 'solvent': 
            solvent_num += 1 
        else: 
            raise Warning('mol_type is not set') 
    return solvent_num 


def count_solute(*args: object) -> int:
    ''' 
    @Description: 
        return the number of solute molecules 
    '''
    solute_num = 0 
    for mol in args:
        if mol.mol_type == 'solute': 
            solute_num += 1 
        else: 
            raise Warning('mol_type is not set') 
    return solute_num