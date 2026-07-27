import sys
import os
import subprocess
import logging
from abc import * 
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
        self.logger.setLevel(logging.INFO)
        self.logger.propagate = False

        formatter = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")

        if not any(isinstance(h, logging.FileHandler) for h in self.logger.handlers):
            file_handler = logging.FileHandler(filename = "autosolvate.log", mode = "a", encoding="utf-8")
            file_handler.setFormatter(formatter)
            file_handler.setLevel(logging.INFO)
            self.logger.addHandler(file_handler)

        # Keep console output minimal: only WARNING+ goes to stdout.
        if not any(isinstance(h, logging.StreamHandler) and not isinstance(h, logging.FileHandler) for h in self.logger.handlers):
            stream_handler = logging.StreamHandler(sys.stdout)
            stream_handler.setFormatter(formatter)
            stream_handler.setLevel(logging.WARNING)
            self.logger.addHandler(stream_handler)

        # Best-effort place to store output file paths produced by this docker.
        # Subclasses can populate this (e.g., in predict_output) to enable concise milestone output.
        self.output_files = []
        os.makedirs(self.workfolder, exist_ok=True)

    @staticmethod
    def resolve_executable(executable: str, amberhome: str = None) -> str:
        """Return executable, preferring an AMBERHOME/bin override when provided."""
        if not amberhome:
            return executable
        expanded = os.path.expandvars(amberhome).rstrip("/")
        if os.path.basename(expanded) == "bin":
            base = expanded
        else:
            base = os.path.join(expanded, "bin")
        return os.path.join(base, executable)

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
        print(f"[AutoSolvate] Running {self.__class__.__name__.replace('Docker', '')}", flush=True)
        print(f"[AutoSolvate] CMD: {cmd}", flush=True)
        if not isinstance(self.exeoutfile, str) or not self.exeoutfile:
            exeout = sys.stdout
            self.logger.info("CMD: {}".format(cmd))
            subprocess.run(cmd, shell = True, stdout=exeout, stderr=sys.stderr)
        else:
            exeout = open(self.exeoutfile, "a")
            self.logger.info("CMD: {}".format(cmd))
            subprocess.run(cmd, shell = True, stdout=exeout, stderr=sys.stderr)
            exeout.close()
        produced = []
        if isinstance(getattr(self, "output_files", None), (list, tuple)):
            produced = [str(p) for p in self.output_files if p]
        if produced:
            existing = [p for p in produced if os.path.exists(p)]
            if existing:
                existing_str = "\n".join(existing)
                print(
                    f"[AutoSolvate] Done  {self.__class__.__name__.replace('Docker', '')} Generated Files:\n{existing_str}",
                    flush=True,
                )
            else:
                print(f"[AutoSolvate] Done  {self.__class__.__name__.replace('Docker', '')}", flush=True)
        else:
            print(f"[AutoSolvate] Done  {self.__class__.__name__.replace('Docker', '')}", flush=True)
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
