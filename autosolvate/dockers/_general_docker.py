import sys
import os
import subprocess
import logging
from abc import * 
print(__file__, __package__)
from ..Common import *

time_execute = 0

class GeneralDocker(ABC):
    """Universal docker template, cannot be instantiated"""

    def __init__(self, 
                 executable:            str = "",
                 workfolder:            str = WORKING_DIR,
                 exeoutfile:            str = None,
    ) -> None:  
        self.executable      = executable
        self.workfolder      = os.path.abspath(workfolder)
        self.exeoutfile      = os.path.abspath(exeoutfile) if isinstance(exeoutfile, str) else exeoutfile
        self.logger                     = logging.getLogger(name = "GeneralDocker")
        self.output_handler                = logging.FileHandler(filename = "autosolvate.log", mode = "a", encoding="utf-8")
        # self.output_handler             = logging.StreamHandler()
        self.output_formater            = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
        self.output_handler.setFormatter(self.output_formater)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(self.output_handler)
        os.makedirs(self.workfolder, exist_ok=True)

    @abstractmethod
    def check_system(self):
        raise NotImplementedError

    @abstractmethod
    def generate_cmd(self):
        raise NotImplementedError

    @abstractmethod
    def generate_input(self):
        raise NotImplementedError
    
    @abstractmethod
    def predict_output(self):
        raise NotImplementedError 


    def execute(self, cmd):
        cwd = os.getcwd()
        os.chdir(self.workfolder)
        if not isinstance(self.exeoutfile, str) or not self.exeoutfile:
            exeout = sys.stdout
            self.logger.info("CMD: {}".format(cmd))
            subprocess.run(cmd, shell = True, stdout=exeout, stderr=sys.stdout)
        else:
            exeout = open(self.exeoutfile, "w")
            self.logger.info("CMD: {}".format(cmd))
            subprocess.run(cmd, shell = True, stdout=exeout, stderr=sys.stdout)
            exeout.close()
        os.chdir(cwd)

    @abstractmethod
    def check_output(self):
        raise NotImplementedError
    
    @abstractmethod
    def process_output(self):
        raise NotImplementedError
    
    @abstractmethod
    def run(self):
        raise NotImplementedError
