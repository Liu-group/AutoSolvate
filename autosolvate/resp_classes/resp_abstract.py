"""
resp.py

Generate resp charge with various quantum chemistry packages

Handles the primary functions
"""
import copy
import logging
import os
import re
import subprocess
from abc import ABC
from typing import Any, Dict, Optional

from openbabel import pybel

from autosolvate.utils.tools import get_residue_name_from_pdb

RESP_ALLOWED_KEYS = {
    "pdbfile",
    "molname",
    "qm_program",
    "qm_exe",
    "qm_dir",
    "rundir",
    "xyzfile",
    "nprocs",
    "nnodes",
    "ncpus",
    "srun_use",
    "charge",
    "spinmult",
    "gamessversion",
    "logger_name",
    "logger",
    "residue_name",
}

AVAILABLE_QM_PROGRAMS = {"gaussian", "gamess", "orca"}
_RESIDUE_NAME_PATTERN = re.compile(r"^[A-Za-z0-9]{1,4}$")


class RespABC(ABC):
    ### Attributes that are in common for different QM packages should go in here.
    def __init__(self, **kwargs):
        parent_logger = kwargs.get("logger")
        parent_logger_name = kwargs.get("logger_name")
        self.logger = self._initialize_logger(parent_logger, parent_logger_name)
        self.logger.debug("RESP inputs: %s", kwargs)
        self.keywords_avail = RESP_ALLOWED_KEYS
        self.pdbfile = False
        self.molname = kwargs.get("molname", "undef")
        self.qm_program = kwargs.get("qm_program", "gamess").lower()
        self.qm_exe = kwargs.get("qm_exe")
        self.qm_dir = kwargs.get("qm_dir")
        self.rundir = kwargs.get("rundir", os.getcwd())
        self.rundir = os.path.abspath(self.rundir)
        os.makedirs(self.rundir, exist_ok=True)
        self.resp_scr_dir = os.path.join(os.path.abspath(self.rundir), 'resp_scr')
        self.initialization_check(**kwargs)
        self.residue_name = self._resolve_residue_name(kwargs.get("residue_name"))
        self.logger.debug("RESP residue name set to %s", self.residue_name)
        self.xyzfile = kwargs.get("xyzfile", "undef")
        self.nprocs = kwargs.get("nprocs")
        if isinstance(self.nprocs, str) and self.nprocs.isdigit():
            self.nprocs = int(self.nprocs)

    def initialization_check(self, **kwargs):
        for key in kwargs:
            if key not in self.keywords_avail:
                self.logger.warning("RESP charge fitting: unrecognized key '%s'", key)

        if self.qm_program not in AVAILABLE_QM_PROGRAMS:
            raise ValueError(f"RESP charge fitting: unsupported program '{self.qm_program}'")

        if "pdbfile" not in kwargs:
            raise KeyError("RESP charge fitting requires 'pdbfile'")
        self.pdbfile = kwargs["pdbfile"]
        self.charge = kwargs.get("charge", 0)
        self.spinmult = kwargs.get("spinmult", 1)
        self.molecule = self.pdb2obmol(self.pdbfile, charge=self.charge, spinmult=self.spinmult)
        if self.molname == "undef":
            self.molname = self.pdb2molname(self.pdbfile)
        self.command_log = os.path.join(self.rundir, f"{self.molname}_{self.qm_program}_commands.log")

    def pdb2obmol(self, pdbfile, charge=0, spinmult=1):
        molecule = pybel.readfile('pdb', pdbfile).__next__()
        obmol = molecule.OBMol
        if charge != 0:
           obmol.SetTotalCharge(charge)
        if spinmult > 1:
           obmol.SetTotalSpinMultiplicity(spinmult)
        return obmol
     
    def pdb2molname(self, pdbfile):
        if pdbfile.find('.pdb') == -1:
             raise ValueError('The input argument must be a file name ended with .pdb')
        molname = copy.copy(pdbfile)
        molname = molname.replace('.pdb','')
        molname = molname.replace('.xyz','')
        return molname

    def run_shell_command(
        self,
        cmd: str,
        description: str = None,
        log_file: str = None,
        produced_file: Optional[str] = None,
        purpose: Optional[str] = None,
    ):
        if description:
            self.logger.info(description)
        self.logger.info("CMD: %s", cmd)
        logfile = log_file or self.command_log
        log_dir = os.path.dirname(logfile)
        if log_dir:
            os.makedirs(log_dir, exist_ok=True)
        with open(logfile, "a", encoding="utf-8") as handle:
            subprocess.check_call(cmd, shell=True, stdout=handle, stderr=subprocess.STDOUT)
        if produced_file and purpose:
            self.logger.info("The %s is generated at %s", purpose, os.path.abspath(produced_file))

    def _initialize_logger(self, parent_logger: Optional[logging.Logger], parent_logger_name: Optional[str]):
        def _clone_handlers(src_logger: logging.Logger) -> logging.Logger:
            child = logging.getLogger(f"{src_logger.name}.{self.__class__.__name__}")
            child.setLevel(src_logger.level)
            child.propagate = False
            for handler in src_logger.handlers:
                if handler not in child.handlers:
                    child.addHandler(handler)
            return child

        if parent_logger and parent_logger.handlers:
            return _clone_handlers(parent_logger)

        if parent_logger_name:
            parent_logger = logging.getLogger(parent_logger_name)
            if parent_logger.handlers:
                return _clone_handlers(parent_logger)

        logger = logging.getLogger(self.__class__.__name__)
        if not logger.handlers:
            handler = logging.FileHandler("autosolvate.log", mode="a", encoding="utf-8")
            formatter = logging.Formatter("%(asctime)s %(name)s %(levelname)s: %(message)s", "%H:%M:%S")
            handler.setFormatter(formatter)
            logger.addHandler(handler)
            logger.setLevel(logging.INFO)
        return logger

    def _resolve_residue_name(self, supplied: Optional[str]) -> str:
        candidate = (supplied or "").strip()
        if candidate:
            return self._validate_residue_name(candidate)
        derived = ""
        if self.pdbfile and os.path.isfile(self.pdbfile):
            try:
                derived = get_residue_name_from_pdb(self.pdbfile) or ""
            except Exception as exc:
                self.logger.warning("Failed to derive residue name from %s: %s", self.pdbfile, exc)
        if derived:
            self.logger.info("Using residue name '%s' derived from PDB", derived)
            return self._validate_residue_name(derived)
        self.logger.warning("Residue name not provided; defaulting to MOL")
        return "MOL"

    def _validate_residue_name(self, residue: str) -> str:
        normalized = residue.strip().upper()
        if not normalized:
            raise ValueError("Residue name is empty after stripping whitespace")
        if not _RESIDUE_NAME_PATTERN.match(normalized):
            raise ValueError("Residue name '%s' is invalid; use 1-4 alphanumeric characters" % normalized)
        return normalized