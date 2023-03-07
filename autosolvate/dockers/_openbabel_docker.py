import getopt, sys, os, subprocess
from typing import Iterable

from ..molecule import *
from ._general_docker import GeneralDocker




class OpenBabelDocker(GeneralDocker):

    def __init__(self, 
                 inp_format:            str = 'pdb',
                 out_format:            str = "mol2",
                 workfolder:            str = WORKING_DIR,
                 exeoutfile:            str = None,
    ) -> None:
        super().__init__(
            executable = '$AMBERHOME/bin/parmchk2',
            workfolder = workfolder,
            exeoutfile = exeoutfile)
        self.inp_format     = inp_format
        self.out_format     = out_format
        self.logger.name    = self.__class__.__name__
        self.outfile        = ""
        self.cmds           = []

    def check_system(self, mol: System) -> None:
        if getattr(mol, self.inp_format, None) is None:
            self.logger.critical("The {} file is not found!".format(self.inp_format))
            raise FileNotFoundError("The {} file is not found!".format(self.inp_format))
        self.logger.info("The {} file is found: {}".format(self.inp_format, getattr(mol, self.inp_format)))
        if isinstance(self.out_format, str):
            self.logger.info("Will generate {} file".format(self.out_format, ))
        elif isinstance(self.out_format, Iterable):
            self.logger.info("Will generate {} file".format(self.out_format, ))

    def predict_output(self, mol: System) -> None:
        if isinstance(self.out_format, str):
            if mol.check_exist(self.out_format):
                self.logger.warning("The {} file {} will be overwritten!".format(self.out_format, getattr(mol, self.out_format)))
            self.logger.info("Will generate {} file at {}".format(self.out_format, getattr(mol, self.out_format)))
            self.outfile = mol.reference_name + "." + self.out_format
        elif isinstance(self.out_format, Iterable):
            self.outfile = []
            for fmt in self.out_format:
                if mol.check_exist(fmt):
                    self.logger.warning("The {} file {} will be overwritten!".format(fmt, getattr(mol, fmt)))
                self.logger.info("Will generate {} file at {}".format(fmt, getattr(mol, fmt)))
                self.outfile.append(mol.reference_name + "." + fmt)

    def generate_input(self, mol: System) -> None:
        pass

    def generate_cmd(self, mol:System) -> str:
        if isinstance(self.out_format, str) and isinstance(self.outfile, str):
            cmd = 'obabel -i {} {} -o {} -O {}'.format(self.inp_format, getattr(mol, self.inp_format), self.out_format, self.outfile)
            self.cmds.append(cmd)
        elif isinstance(self.out_format, Iterable) and isinstance(self.outfile, Iterable):
            for fmt, outfile in zip(self.out_format, self.outfile):
                cmd = 'obabel -i {} {} -o {} -O {}'.format(self.inp_format, getattr(mol, self.inp_format), fmt, outfile)
                self.cmds.append(cmd)
        else:
            self.logger.critical("The output format and output file are not consistent!")
            raise RuntimeError("The output format and output file are not consistent!")
        return cmd

    def execute(self, cmd:str) -> None:
        self.logger.info("Running OpenBabel ...")
        cwd = os.getcwd()
        os.chdir(self.workfolder)
        if not isinstance(self.exeoutfile, str) or not self.exeoutfile:
            exeout = sys.stdout
            for cmd in self.cmds:
                self.logger.info("CMD: {}".format(cmd))
                subprocess.run(cmd, shell = True, stdout=exeout, stderr=sys.stdout)
        else:
            exeout = open(self.exeoutfile, "w")
            for cmd in self.cmds:
                self.logger.info("CMD: {}".format(cmd))
                subprocess.run(cmd, shell = True, stdout=exeout, stderr=sys.stdout)
            exeout.close()
        os.chdir(cwd)
        # obConversion = ob.OBConversion() 
        # obConversion.SetInAndOutFormats("xyz", "pdb") 
        # OBMOL = ob.OBMol() 
        # obConversion.ReadFile(OBMOL, mol.xyz) 
        # obConversion.WriteFile(OBMOL, mol.name+'/'+mol.name+'.pdb')

    def check_output(self):
        if isinstance(self.out_format, str) and isinstance(self.outfile, str):
            if not os.path.exists(self.outfile):
                self.logger.critical("The {} file {} is not found!".format(self.out_format, self.outfile))
                raise FileNotFoundError("The {} file {} is not found!".format(self.out_format, self.outfile))
            self.logger.info("OpenBabel successfully generated the {} file: {}".format(self.out_format, self.outfile))
        elif isinstance(self.out_format, Iterable) and isinstance(self.outfile, Iterable):
            for fmt, outfile in zip(self.out_format, self.outfile):
                if not os.path.exists(outfile):
                    self.logger.critical("The {} file {} is not found!".format(fmt, outfile))
                    raise FileNotFoundError("The {} file {} is not found!".format(fmt, outfile))
                self.logger.info("OpenBabel successfully generated the {} file: {}".format(fmt, outfile))
    
    def process_output(self, mol: System) -> None:
        if isinstance(self.out_format, str) and isinstance(self.outfile, str):
            setattr(mol, self.out_format, self.outfile)
        elif isinstance(self.out_format, Iterable) and isinstance(self.outfile, Iterable):
            for fmt, outfile in zip(self.out_format, self.outfile):
                setattr(mol, fmt, outfile)
        mol.update()

    def run(self, mol: System) -> None:
        self.logger.name = self.__class__.__name__
        self.check_system(mol)
        self.predict_output(mol)
        self.generate_input(mol)
        cmd = self.generate_cmd(mol)
        self.execute(cmd)
        self.check_output()
        self.process_output(mol)

