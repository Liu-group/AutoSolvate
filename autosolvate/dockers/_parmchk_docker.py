import getopt, sys, os, subprocess

from ..Common import * 
from ..molecule import *
from ..utils import srun
from ._general_docker import GeneralDocker


class ParmchkDocker(GeneralDocker):
    def __init__(self, 
                 out_format:            str = 'frcmod',
                 workfolder:            str = WORKING_DIR,
                 exeoutfile:            str = None,
    ) -> None:
        super().__init__(
            executable = '$AMBERHOME/bin/parmchk2',
            workfolder = workfolder,
            exeoutfile = exeoutfile)
        self.out_format     = out_format
        self.logger.name    = self.__class__.__name__
        self.outfile        = ""

    def get_output_name(self, mol:System, fmt:str) -> str:
        frcmod = ".".join([mol.reference_name, self.out_format])
        return frcmod

    def set_executable(self) -> str:
        return '$AMBERHOME/bin/parmchk2'

    def set_input(self, mol: System) -> str:    
        if mol.mol2 is not None:
            return '-i %s -f %s' % (mol.mol2, 'mol2') 
        raise Exception('only support mol2 as input format and {} does not exist'.format(mol.mol2)) 

    def set_output(self, mol: object) -> str: 
        return '-o %s' % (self.outfile)

    def check_system(self, mol:System):
        self.logger.info("Checking system {:s}...".format(mol.name))
        if mol.check_exist("mol2"):
            self.logger.info("System mol2: {:s}".format(mol.mol2))
        else:
            self.logger.critical("Only support mol2 as input format but {} does not exist.".format(mol.mol2))
            raise RuntimeError("Only support mol2 as input format but {} does not exist.".format(mol.mol2))
    
    def generate_input(self, mol:System):
        return 
    @srun()
    def generate_cmd(self, mol:System):
        '''
        @EXAMPLE: 
        $AMBERHOME/bin/parmchk2 -i 1.prmtop -f mol2 -o 1.frcmod
        '''
        cmd =  self.executable          + ' ' 
        cmd += self.set_input(mol)      + ' ' 
        cmd += self.set_output(mol)     + ' ' 
        return cmd 
    
    def predict_output(self, mol:System):
        outname = self.get_output_name(mol, self.out_format)
        self.logger.info("The {} file will be generated at {}".format(self.out_format, outname))
        self.outfile = outname
        if os.path.exists(outname):
            self.logger.warn("Found a existing file with the same name: {}".format(outname))
            self.logger.warn("This file will be Overwritten!".format(outname))
        return outname
    
    def check_output(self, mol: System):
        success = True
        outfile = self.outfile
        if os.path.exists(outfile):
            self.logger.info("Successfully generated target {} file: {}".format(self.out_format, outfile))
            success = True
        else:
            self.logger.critical("Failed to generate target {} file: {}".format(self.out_format, outfile))
            success = False
        if not success:
            raise RuntimeError("Failed to generate target {} file: {}".format(self.out_format, outfile))
    
    def process_output(self, mol:System):
        frcmodfile = self.outfile
        setattr(mol, self.out_format, frcmodfile)
        mol.update()

    def run(self, mol:System):
        self.logger.name = self.__class__.__name__
        self.check_system(mol)
        self.predict_output(mol)
        self.generate_input(mol)
        cmd = self.generate_cmd(mol)
        self.execute(cmd)
        self.check_output(mol)
        self.process_output(mol)

