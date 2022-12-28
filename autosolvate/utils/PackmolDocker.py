from openbabel import pybel
from openbabel import openbabel as ob
from multipledispatch import dispatch 
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools

class PackmolDocker:

    def __init__(self): 
        '''
        @NOTE: 
        1. I try to not let a function oriented class initialize any variables 
        2. Do not initialize a list as instance variable, because it will be shared by all instances 
           (python will only create one list for all instances)
        '''
        self.out_name                   = 'system.pdb'   


    def run(self, box: object) -> None:
        check_mol(*box.solute_list)
        check_mol(*box.solvent_list) 
        r'''
        @TODO 
        1. finish this function with parameter box  
        '''
        self.write_packmol_inp(box)
        cmd = self.generate_cmd()
        tools.submit(cmd)
        

    @tools.srun()
    def generate_cmd(self) -> str:
        cmd = 'packmol < packmol.inp > packmol.log'
        return cmd 
    

    def write_packmol_inp(self, box: object) -> None: 
        r'''
        @TODO 
        1. finish this function 
        '''
        f = open('packmol.inp', 'w') 
        f.write('{:<15} {:<3.2f}                \n'.format('tolerance', box.closeness))
        f.write('{:<15} {:<5}                   \n'.format('filetype', 'pdb'))
        f.write('{:<15} {:<5}                   \n'.format('output', self.out_name)) 
        self.load_solute(f,  box)
        self.load_solvent(f, box)
        f.close()


    def load_solute(self, doc: object, box: object) -> None: 
        r'''
        @TODO
        1. support multiple solutes
        2. can be simplified further 

        @QUESTION
        1. does have to use 'solute.pdb' here? 
        2. check if {com} is correct
        '''
        if len(box.solute_list) == 0:
            raise Exception('no solute found')
        
        if len(box.solute_list) > 1:
            raise Exception('multiple different solutes not supported yet')
        
        if len(box.solute_list) == 1:
            solute      = solute_list[0]
            solute_pos  = box.cubesize / 2.0 

            doc.write('{}           \n'.format('# add the solute'))
            doc.write('{:<15} {:<5} \n'.format('structure', solute.pdb))
            doc.write('{:<15} {:<5} \n'.format('number',    box.duplicate_solute_num))

            if box.duplicate_solute_num == 1:
                doc.write('{:<10} {pos} {pos} {pos} {com} {com} {com}\n'.format('fixed', pos=solute_pos, com='0.'))
                doc.write('{:<15} {:<5}                              \n'.format('resnumbers', '2'))
                doc.write('{:<15}                                    \n'.format('centerofmass'))
                doc.write('{:<15}                                    \n'.format('end structure'))

            else:
                doc.write('{:<10} {pos} {pos} {pos} {cube}           \n'.format('inside cube', pos=solute_pos, cube=box.cubesize))
                doc.write('{:<15} {:<5}                              \n'.format('resnumbers', '2'))
            doc.write('{:<15}                                        \n'.format('end structure'))
            doc.write('\n')
        


    def load_solvent(self, doc: object, box: object) -> None: 
        
        if len(box.solvent_list) == 0: 
            raise Exception('no solvent found')
        
        for solvent in box.solvent_list:
            doc.write('{}                                       \n'.format('# add the solvent'))
            doc.write('{:<15} {:<5}                             \n'.format('structure', solvent.pdb))
            doc.write('{:<15} {:<5}                             \n'.format('number',    box.duplicate_solvent_num))
            doc.write('{:<10} {pos} {pos} {pos} {boxsize}       \n'.format('inside cube', pos=box.cubesize/2.0, boxsize=box.cubesize))
            doc.write('{:<15}                                   \n'.format('end structure'))
            doc.write('\n')    
        


#METHODS 
def check_mol(*args: object) -> None:
    for mol in args: 
        if mol.pdb is None: 
            raise Exception('pdb file is not loaded')
        if mol.mol_type is None: 
            raise Exception('mol_type is not loaded') 


