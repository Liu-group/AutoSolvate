from openbabel import pybel
from openbabel import openbabel as ob
from Common import * 
import numpy as np
import getopt, sys, os, subprocess
import tools

class AntechamberWapper():

    SUPPORT_INPUT_FORMATS           = ['pdb']
    SUPPORT_CHARGE_FITTING_METHODS  = ['bcc']

    def __init__(self, inputfile, charge=0, multiplicity=1, charge_method='bcc', output_format='mol2', residue_name='MOL'):      
        r'''
        @EXAMPLE: 
        >>> AntechamberWapper('test.pdb', charge=0, multiplicity=1, charge_method='bcc', output_format='mol2', residue_name='MOL').run() 
        $AMBERHOME/bin/antechamber -i test.pdb -fi pdb -o test.mol2 -fo mol2   -c bcc  

        
        >>> AntechamberWapper('test.pdb', charge=1, multiplicity=2, charge_method='bcc', output_format='mol2', residue_name='SLU').run()
        $AMBERHOME/bin/antechamber -i test.pdb -fi pdb -o test.mol2 -fo mol2 -nc 1 -m 2 -c bcc -rn SLU 
        '''
        #input
        self.inputfile              = inputfile
        #output
        self.output_format          = output_format

        self.charge                 = charge
        self.multiplicity           = multiplicity 
        self.charge_method          = charge_method 
        self.residue_name           = residue_name 

    def set_executable(self):
        return '$AMBERHOME/bin/antechamber'
        

    def set_input(self):         
        basename, name, ext = tools.get_basename_name_ext(self.inputfile) 
        
        if ext not in self.SUPPORT_INPUT_FORMATS: 
            raise Exception('Input file format not supported in class AntechamberWapper()') 
        return '-i %s -fi %s' % (self.inputfile, ext) 


    def set_output(self):
        basename, name, ext = tools.get_basename_name_ext(self.inputfile) 
        return '-o %s.mol2 -fo mol2' % name


    def set_charge(self):
        if self.charge == 0:
            return '' 
        return '-nc %d' % self.charge 


    def set_multiplicity(self): 
        if self.multiplicity == 1:
            return '' 
        return '-m %d' % self.multiplicity 
        

    def set_charge_method(self): 
        if self.charge_method not in self.SUPPORT_CHARGE_FITTING_METHODS:
            raise Exception('Charge fitting method not supported in class AntechamberWapper()')
        return '-c %s' % self.charge_method


    def set_residue_name(self): 
        if self.residue_name == 'MOL':
            return '' 
        return '-rn %s' % self.residue_name

    @tools.srun()
    def generate_cmd(self):
        cmd =  self.set_executable()        + ' '
        cmd += self.set_input()             + ' '
        cmd += self.set_output()            + ' ' 
        cmd += self.set_charge()            + ' ' 
        cmd += self.set_multiplicity()      + ' ' 
        cmd += self.set_charge_method()     + ' ' 
        cmd += self.set_residue_name()      + ' ' 
        return cmd 

    def run(self):
        cmd = self.generate_cmd()
        if DRY_RUN:
            print(cmd) 
            return
        else:
            subprocess.call(cmd, shell=True) 





if __name__ == '__main__': 
    import doctest

    global DRY_RUN 
    DRY_RUN = True

    doctest.testmod() 
    