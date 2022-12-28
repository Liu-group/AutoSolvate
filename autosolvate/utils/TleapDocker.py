from openbabel import pybel
from openbabel import openbabel as ob
from multipledispatch import dispatch 
import numpy as np
import getopt, sys, os, subprocess
from Common     import *
from Molecule   import Molecule 
from SolventBox import SolventBox  
import tools


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
        self._load_ions             = False
        

    def run(self, o: object) -> None:

        if isinstance(o, Molecule):
            mol = o 
            check_mol_arributes(mol)
            self.write_tleap_in(mol)
            cmd = self.generate_cmd()
            tools.submit(cmd)
            return 
        
        if isinstance(o, SolventBox):
            '''
            @TODO
            modify 
            finish this function
            '''
            box = o 
            check_mol_arributes(*box.solute_list)
            check_mol_arributes(*box.solvent_list)

            if len(box.solute_list) == 1 and len(box.solvent_list) == 1:
                '''
                @NOTE: 
                1. this is for one solute and one solvent 
                '''
                solute  = box.solute_list[0]
                solvent = box.solvent_list[0]
                self.write_tleap_in(solute, solvent, box.closeness, box.cubesize)
                cmd = self.generate_cmd()
                tools.submit(cmd)
            else: 
                raise Exception('now only support one solute and one solvent') 
            return

        raise Warning('invalid input type')
        

    @tools.srun()
    def generate_cmd(self) -> str:
        cmd = 'tleap -s -f leap.in > leap.out'
        return cmd


    #dispatched write_tleap_in function 
    @dispatch(object)
    def write_tleap_in(self, mol: object) -> None:  
        r'''
        @QUESTION 
        1. Is check mol necessary?, what is it for? 
        2. Does self count as an object? in @dispatch(object) 

        @TODO 
        1. check if 20 pt is enough for all strings 
        '''
        f = open('leap.in', 'w')
        self.load_forcefield(f)
        self.load_mol(f, mol, frcmod=True, mol2=True)  
        f.write('{:<20} {:<5} {:<15}        \n'.format('saveoff', mol.residue_name, mol.name + '.lib'))
        f.write('{:<20} {:<5} {:<15}        \n'.format('savepdb', mol.residue_name, mol.name + '.pdb'))
        f.write('{:<20} {:<5} {:<15} {:<15} \n'.format('saveamberparm', mol.residue_name, mol.name + '.prmtop', mol.name + '.inpcrd')) 
        f.write('{:<20}                     \n'.format('quit'))
        f.close() 

    
    @dispatch(object, object, float, int)
    def write_tleap_in(self, 
                       solute:      object, 
                       solvent:     object,
                       closeness:   float, 
                       cubsize:     int 
    ) -> None: 
        r'''
        @TODO: 
        1. combine change pass in closeness and cubsize to pass in solventbox object 

        @QUESTION: 
        1. does solvent need mol2 file? 
        '''
        #setting
        if solute.charge != 0:
            self._load_ions = True 

        #set pos 
        pos = cubsize / 2.0 

        #output name 
        output_name = solute.name + '_solvated'

        #write tleap.in 
        f = open('leap.in', 'w')
        self.load_forcefield(f) 
        self.load_mol(f, solvent, frcmod=True, lib=True)
        self.load_mol(f, solute,  frcmod=True, mol2=True, check_mol=True)
        self.load_solventbox(f, solute, solvent, pos, closeness)
        if self._load_ions:
            self.load_ions(f, solute, solvent)
        f.write('{:<20}  {:<20}         \n'.format('savepdb mol', output_name+'.pdb')) 
        f.write('{:<20}  {:<20}  {:<20} \n'.format('saveamberparm mol', output_name+'.prmtop', output_name+'.inpcrd'))
        f.write('{:<20}                 \n'.format('quit'))      
    #end of write_tleap_in function 


    def load_forcefield(self, doc: object) -> None:
        doc.write('{:<20}  {:<20}       \n'.format('source', self._ff14SB)) 
        doc.write('{:<20}  {:<20}       \n'.format('source', self._gaff))
        if self._load_ions:
            doc.write('{:<20}  {:<20}   \n'.format('source', self._water_tip3p))
        

    def load_mol(self,  
                 doc:           object, 
                 mol:           object, 
                 frcmod:        bool = False, 
                 mol2:          bool = False, 
                 lib:           bool = False, 
                 check_mol:     bool = False
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


    def load_solventbox(self, 
                        doc:            object, 
                        solvent:        object, 
                        pos:            float, 
                        closeness:      float
    ) -> None:
        if solvent in AMBER_SOLVENT_LIST:
            box = solvent.box
            
        else:
            box = solvent.name 

        doc.write('{:<20}  {:<10}  {:<10}  {:<3}  {:<5} \n'.format(
            'solvatebox mol', box, pos, 'iso', closeness
        ))
    

    def load_ion(self, doc: object, solute: object) -> None:
        if solute.charge > 0: 
            ion = 'Cl-'
        elif solute.charge < 0: 
            ion = 'Na+' 
        else: 
            raise Exception('invalid charge') 
        doc.write('{:<20} {:<5} {:<5} \n'.format('addIons2 mol', ion + ''))
        doc.write('{:<20}             \n'.format('check mol')) 


    def write_packmol(self, closeness: float, *args: object) -> None: 
        '''
        @TODO: 
            1. support multiple solute and solvent in the future             
        '''


#METHODS 
def check_mol_arributes(*args: object) -> None: 
    '''
    @Description: 
        check if important attributes of the molecule are set for tleap 
    ''' 
    for mol in args:
        if mol.mole_type is None: 
            raise Exception('mole_type is not set') 
        if mol.residue_name is None:
            raise Exception('residue name is not set')     
        if mol.mol2 is None:
            raise Exception('mol2 file is not set')
        if mol.frcmod is None:
            raise Exception('frcmod file is not set')
