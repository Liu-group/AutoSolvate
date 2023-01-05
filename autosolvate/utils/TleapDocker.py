from openbabel import pybel
from openbabel import openbabel as ob
from multipledispatch import dispatch 
import numpy as np
import getopt, sys, os, subprocess
from Common     import *
from Molecule   import Molecule, AMBER_SOLVENT_LIST
from SolventBox import SolventBox  
import tools


class TleapDocker: 

    _ff14SB         = 'leaprc.protein.ff14SB'       #Source leaprc file for ff14SB protein force field
    _gaff           = 'leaprc.gaff'                 #Source leaprc file for gaff force field
    _water_tip3p    = 'leaprc.water.tip3p'


    def __init__(self):     
        '''
        @NOTE:
        1. do not store any molecule obejct in this class attribute 
           because this class is a function class, not a data class 
        
        @TODO: 
        1. check if the logic of _load_ions is correctly applied
        '''
        self._load_ions             = False


    def run(self, o: object) -> None:
        
        #Molecule as input 
        if isinstance(o, Molecule):
            os.chdir(o.name)
            mol = o 
            check_mol(mol)
            self.write_tleap_in(mol)
            cmd = self.generate_cmd()
            tools.submit(cmd)
            os.chdir(WORKING_DIR)
            return 
        
        
        #SolventBox as input 
        if isinstance(o, SolventBox):
            os.chdir(o.name)
            '''
            @TODO
            1. logic for checking these conditions can be simplified or better written 
            2. need to discuss with Dr Liu about the logic of this function 
            '''
            box = o 
            check_mol(*box.solute_list)
            check_mol(*box.solvent_list)

            if len(box.solute_list) == 1 and len(box.solvent_list) == 1:
                '''
                @NOTE: 
                1. this is for one solute and one solvent 
                '''
                solute      = box.solute_list[0]
                solvent     = box.solvent_list[0]
                system_pdb  = box.system_pdb
                
                if system_pdb is not None:
                    self.write_tleap_in(solute, solvent, box.closeness, box.cubesize, box.system_pdb)
                
                elif solvent in AMBER_SOLVENT_LIST: 
                    self.write_tleap_in(solute, solvent, box.closeness, box.cubesize)
                
                elif solvent.frcmod is not None and solvent.lib is not None:
                    self.write_tleap_in(solute, solvent, box.closeness, box.cubesize)

                else: 
                    raise Exception('case not considered yet') 

                cmd = self.generate_cmd()
                tools.submit(cmd)
            else: 
                raise Exception('only support one solute and one solvent') 
            
            os.chdir(WORKING_DIR)
            return

        raise Warning('invalid input type')
        

    @tools.srun()
    def generate_cmd(self) -> str:
        cmd = 'tleap -s -f leap.in > leap.out'
        return cmd










    #dispatched write_tleap_in function 
    #CASE 1. only one mol 
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

    #CASE 2. one solute and one solvent 
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

        pos         = cubsize / 2.0 
        output_name = solute.name + '_solvated'

        #write tleap.in 
        f = open('leap.in', 'w')

        self.load_forcefield(f)
        self.load_mol(f, solvent, frcmod=True, lib=True)
        self.load_mol(f, solute,  frcmod=True, mol2=True, check=True)
        self.load_solventbox(f, solute, solvent, pos, closeness)
        if self._load_ions:
            self.load_ions(f, solute, solvent)
        
        f.write('{:<20}  {:<5}   {:<15}         \n'.format('savepdb', solute.residue_name, output_name+'.pdb')) 
        f.write('{:<20}  {:<5}   {:<20}  {:<20} \n'.format('saveamberparm', solute.residue_name, output_name+'.prmtop', output_name+'.inpcrd'))
        f.write('{:<20}                         \n'.format('quit'))      
        f.close()
    
    
    #CASE 3. one solute, one solvent, and  system_pdb
    @dispatch(object, object, float, int, str) 
    def write_tleap_in(self, 
                       solute:      object, 
                       solvent:     object,
                       closeness:   float, 
                       cubsize:     int, 
                       system_pdb:  str 
    ) -> None: 

        #setting
        if solute.charge != 0:
            self._load_ions = True 

        pos         = cubsize / 2.0 
        output_name = solute.name + '_solvated'

        #write tleap.in 
        f = open('leap.in', 'w')
        
        self.load_forcefield(f) 
        self.load_mol(f, solvent, frcmod=True, mol2=True, check=True)
        self.load_mol(f, solute,  frcmod=True, mol2=True, check=True)
        '''
        @QUESTION: 
        I dont know how to load a ions in this case 
        '''
        f.write('{:<5}  {:<10}   {:<20}            \n'.format('SYS =', 'loadpdb', system_pdb)) 
        f.write('{:<5}  {:<5}    {:<5}  {:<20}     \n'.format('set', 'SYS', 'box', '{'+str(cubsize)+' '+str(cubsize)+' '+str(cubsize)+'}'))
        f.write('{:<20} {:<5}    {:<20} {:<20}     \n'.format('saveamberparm', 'SYS', output_name+'.prmtop', output_name+'.inpcrd'))
        f.write('{:<20} {:<5}    {:<20}            \n'.format('savepdb', 'SYS', output_name+'.pdb'))
        f.write('{:<5}                             \n'.format('quit'))     
        f.close()    
    #end of write_tleap_in function 
        





    #helper functions that write tleap.in 
    def load_forcefield(self, doc: object) -> None:
        doc.write('{:<20}  {:<20}       \n'.format('source', self._ff14SB)) 
        doc.write('{:<20}  {:<20}       \n'.format('source', self._gaff))
        doc.write('{:<20}  {:<20}       \n'.format('source', self._water_tip3p))
        

    def load_mol(self,  
                 doc:           object, 
                 mol:           object, 
                 frcmod:        bool = False, 
                 mol2:          bool = False, 
                 lib:           bool = False, 
                 check:         bool = False
    ) -> None: 
        r'''
        @QUESTION:
        1. what is check mol for? 
        '''
        if mol in AMBER_SOLVENT_LIST: 
            return 
        if mol2:  
            doc.write('{:<5} {:<15} {:<20}  \n'.format(mol.residue_name, '= loadmol2', mol.mol2))
        if frcmod: 
            doc.write('{:<20}  {:<20}       \n'.format('loadamberparams', mol.frcmod))       
        if lib: 
            doc.write('{:<20}  {:<20}       \n'.format('loadoff', mol.lib)) 
        if check: 
            doc.write('{:<20}  {:<20}       \n'.format('check', mol.residue_name))


    def load_solventbox(self, 
                        doc:            object,
                        solute:         object, 
                        solvent:        object, 
                        pos:            float, 
                        closeness:      float
    ) -> None:
        if solvent in AMBER_SOLVENT_LIST:
            solvent_model= solvent.box
            
        else:
            solvent_model = solvent.name 

        doc.write('{:<20}  {:<10}  {:<10}  {:<10}  {:<3}  {:<5} \n'.format(
            'solvatebox', solute.residue_name, solvent_model, pos, 'iso', closeness
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





#METHODS 
def check_mol(*args: object) -> None: 
    '''
    @Description: 
        check if important attributes of the molecule are set for tleap 
    '''
    for mol in args:
        if mol in AMBER_SOLVENT_LIST:
            continue
        if mol.mol_type is None: 
            raise Exception('{} mole_type is not set'.format(mol.name))
        if mol.residue_name is None:
            raise Exception('{} residue name is not set'.format(mol.name))  
        if mol.frcmod is None:
            raise Exception('{} frcmod file is not set'.format(mol.name)) 
