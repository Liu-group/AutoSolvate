from openbabel import pybel
from openbabel import openbabel as ob
from multipledispatch import dispatch 
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools
from Molecule import Molecule 

class PackmolDocker:

    def __init__(self): 
        '''
        @NOTE: 
        1. I try to not let a function oriented class initialize any variables 
        2. Do not initialize a list as instance variable, because it will be shared by all instances 
           (python will only create one list for all instances)
        '''
        self._solute_num                = 0
        self._solvent_num               = 0
       
        self._solute_list               = None
        self._solvent_list              = None
       
        self._duplicate_solute_num      = 1
        self._duplicate_solvent_num     = 100


    def run(self, *mol: object) -> None:
        r'''
        @TODO 
        1. finish this function 
        '''
        self.check_args(*mol)
        self.write_packmol_inp(self._solute_list, self._solvent_list) 
        cmd = self.generate_cmd()




    @tools.srun()
    def generate_cmd(self) -> str:
        cmd = 'packmol < packmol.inp > packmol.log'
        return cmd 


    def check_args(self, *args: object) -> None: 
        check_mol_arributes(*args)
        self._solute_num        = tools.count_solute(*args)
        self._solvent_num       = tools.count_solvent(*args)
        self._solute_list       = tools.get_list_mol_type(*args, mol_type='solute') 
        self._solvent_list      = tools.get_list_mol_type(*args, mol_type='solvent') 
    

    def write_packmol_inp(self, solute_list, solvent_list, closeness: float = 2.0, cube_size: int = 54) -> None:
        r'''
        @TODO 
        1. finish this function 
        '''
        output_name = 'system.pdb'
        f = open('packmol.inp', 'w') 
        f.write('{:<15} {:<3.2f}    \n'.format('tolerance', closeness))
        f.write('{:<15} {:<5}       \n'.format('filetype', 'pdb'))
        f.write('{:<15} {:<5}       \n'.format('output', output_name))
        self.load_solute(f,  solute_list,  cube_size) 
        self.load_solvent(f, solvent_list, cube_size) 
        f.close()


    def load_solute(self, doc: object, solute_list: object, cube_size: int) -> None:
        r'''
        @TODO
        1. support multiple solutes
        2. can be simplified further 

        @QUESTION
        1. does have to use 'solute.pdb' here? 
        2. check if {com} is correct
        '''
        if len(solute_list) == 0: 
            raise Exception('no solute found')
        
        if len(solute_list) > 1: 
            raise Exception('multiple different solutes not supported yet')
        
        if len(solute_list) == 1:
            solute      = solute_list[0]
            solute_pos  = cube_size / 2.0 

            doc.write('{}           \n'.format('# add the solute'))
            doc.write('{:<15} {:<5} \n'.format('structure', solute.pdb))
            doc.write('{:<15} {:<5} \n'.format('number', self._duplicate_solute_num))
            if self._duplicate_solute_num == 1:
                doc.write('{:<10} {pos} {pos} {pos} {com} {com} {com}\n'.format('fixed', pos=solute_pos, com='0.'))
                doc.write('{:<15} {:<5}                              \n'.format('resnumbers', '2'))
                doc.write('{:<15}                                    \n'.format('centerofmass'))
                doc.write('{:<15}                                    \n'.format('end structure'))

            else:
                doc.write('{:<10} {pos} {pos} {pos} {boxsize}       \n'.format('inside cube', pos=solute_pos, boxsize=cube_size))
                doc.write('{:<15} {:<5}                             \n'.format('resnumbers', '2'))
            doc.write('{:<15}       \n'.format('end structure'))
            doc.write('\n')
        


    def load_solvent(self, doc: object, solvent_list: object, cube_size: int) -> None: 
        if len(solvent_list) == 0: 
            raise Exception('no solvent found')
        for solvent in solvent_list: 
            doc.write('{}           \n'.format('# add the solvent'))
            doc.write('{:<15} {:<5} \n'.format('structure', solvent.pdb))
            doc.write('{:<15} {:<5} \n'.format('number', self._duplicate_solvent_num))
            doc.write('{:<10} {pos} {pos} {pos} {boxsize}       \n'.format('inside cube', pos=cube_size/2.0, boxsize=cube_size))
            doc.write('{:<15}       \n'.format('end structure'))
            doc.write('\n')    
        


#METHODS 
def check_mol_arributes(*args: object) -> None:
    for mol in args: 
        if mol.pdb is None: 
            raise Exception('pdb file is not loaded')
        if mol.mol_type is None: 
            raise Exception('mol_type is not loaded') 


