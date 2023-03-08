import getopt
import sys
import os
import shutil
import subprocess
import pkg_resources
import glob
import logging

from typing import Any, List, Tuple
from dataclasses import dataclass, field, asdict

from ..Common import *
from ..utils import *

from logging import DEBUG, INFO, WARN, WARNING, CRITICAL
logging.basicConfig(level = INFO, force = True, handlers=[])

logger              = logging.getLogger(name = "Molecule")
# output_handler      = logging.FileHandler(filename = "log.txt", mode = "a", encoding="utf-8")
output_handler      = logging.StreamHandler()
output_formater     = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
output_handler.setFormatter(output_formater)
if len(logger.handlers) == 0:
    logger.addHandler(output_handler)

OBABEL = "obabel"

@dataclass
class System(object):
    """
    A most basic class that does not contain molecular structures. 
    Only includes the function of adding, deleting, modifying and checking the corresponding files.
    """
    # constants
    
    # required arguments
    name:           str = field(default=None, init = True)
    # positional arguments
    _FILEATTR_FILEPATH_DICT:dict = field(default_factory=dict, init = False)

    # general structural files
    smi:            str = field(default=None, init = False)
    xyz:            str = field(default=None, init = False)
    pdb:            str = field(default=None, init = False)
    cif:            str = field(default=None, init = False)
    sdf:            str = field(default=None, init = False)
    mol2:           str = field(default=None, init = False)

    # Amber files
    lib:            str = field(default=None, init = False)
    off:            str = field(default=None, init = False)
    prep:           str = field(default=None, init = False)
    frcmod:         str = field(default=None, init = False)
    prmtop:         str = field(default=None, init = False)
    inpcrd:         str = field(default=None, init = False)

    # Gromacs files
    top:            str = field(default=None, init = False)
    itp:            str = field(default=None, init = False)
    gro:            str = field(default=None, init = False)

    def __post_init__(self) -> None:
        self._FILEATTR_FILEPATH_DICT = asdict(self)
        self._FILEATTR_FILEPATH_DICT.pop("_FILEATTR_FILEPATH_DICT")
        self._FILEATTR_FILEPATH_DICT.pop("name")
        # self._FILEATTR_FILEPATH_DICT.pop("folder")
        __dict = System.__dict__
        rmkeys = []
        for key in self._FILEATTR_FILEPATH_DICT:
            if key not in __dict:
                rmkeys.append(key)
        for key in rmkeys:
            self._FILEATTR_FILEPATH_DICT.pop(key)
        if not hasattr(self, "folder"):
            self.folder = WORKING_DIR
        self._init_folder()

    def __setattr__(self, __name: str, __value: Any) -> None:
        if not hasattr(self, "_FILEATTR_FILEPATH_DICT"):
            super(System, self).__setattr__(__name, __value)
            return 
        if isinstance(__value, str):
            if os.path.isfile(__value) and len(__name) <= 6:
                __value = os.path.abspath(__value)
                self._FILEATTR_FILEPATH_DICT[__name] = __value
                logger.info("set the '{:s}' file of system '{:s}' to '{:s}'".format(__name, self.name, __value))
                if not self.check_extension(__name, __value):
                    logger.warn("The '{:s}' file of system '{:s}' does not have the correct extension".format(__name, self.name, __value))
            else:
                if __name in self._FILEATTR_FILEPATH_DICT:
                    self._FILEATTR_FILEPATH_DICT[__name] = __value
                    logger.warn("The '{:s}' file of system '{:s}' is set to a non-existent path '{:s}'".format(__name, self.name, __value))
                else:
                    logger.info("set the '{:s}' attribute of system '{:s}' as {:s}".format(__name, self.name, __value))
        else:
            if __name in self._FILEATTR_FILEPATH_DICT:
                logger.warn("set the '{:s}' file of system '{:s}' as a {:s}".format(__name, self.name, type(__value).__name__))
            else:
                logger.info("set the '{:s}' attribute of system '{:s}' as a {:s}".format(__name, self.name, type(__value).__name__))
        super(System, self).__setattr__(__name, __value)

    def check_extension(self, __name, __value:str) -> bool:
        if not isinstance(__value, str):
            return False
        ext = os.path.splitext(__value)[-1][1:]
        return ext == __name

    def check_exist(self, extensions:List[str]) -> bool:
        flag = 1
        if isinstance(extensions, str):
            extensions = [extensions, ]
        for ext in extensions:
            if ext.startswith(r"."):
                ext = ext[1:]
            if not hasattr(self, ext):
                flag *= False
                logger.debug("The '{:s}' file of system '{:s}' is not defined".format(ext, self.name))
            elif not isinstance(self.__getattribute__(ext), str):
                flag *= False
                logger.debug("The '{:s}' attribute of system '{:s}' does noe exist".format(ext, self.name))
            elif not os.path.isfile(self.__getattribute__(ext)):
                flag *= False
                logger.debug("The '{:s}' file of system '{:s}' does not exist".format(ext, self.name))
            else:
                logger.debug("The '{:s}' file of system '{:s}' is '{:s}'".format(ext, self.name, self.__getattribute__(ext)))
        return bool(flag)
    
    def update(self, overwrite = True) -> None:
        if not os.path.exists(self.folder):
            os.makedirs(self.folder)
        for ext in self._FILEATTR_FILEPATH_DICT:
            if not self.check_exist(ext):
                continue
            srcpath = self.__getattribute__(ext)
            newpath = os.path.join(
                self.folder,
                os.path.basename(srcpath),
            )
            if os.path.exists(newpath):
                if os.path.samefile(srcpath, newpath):
                    logger.debug("The '{:s}' file of system '{:s}' already exists".format(ext, self.name))
                elif overwrite:
                    shutil.copy(srcpath, self.folder)
                    logger.info("Copied '{:s}' file of system '{:s}' to '{:s}'".format(ext, srcpath, newpath))
                    logger.info("File {:s} was rewritten".format(newpath))
                    self.__setattr__(ext, newpath)
                else:
                    logger.info("The '{:s}' file of system '{:s}' already exists at {:s} ".format(ext, self.__getattribute__(ext), newpath))
                    self.__setattr__(ext, newpath)
            else:
                shutil.copy(srcpath, self.folder)
                logger.info("Copied '{:s}' file of system '{:s}' to '{:s}'".format(ext, srcpath, newpath))
                self.__setattr__(ext, newpath)

    def copy2folder(self, overwrite = True):
        self.update(overwrite=overwrite)

    def _init_folder(self) -> None:
        if self.folder == WORKING_DIR: 
            logger.info("Use the current working directory {:s}".format(self.folder))
        elif os.path.exists(self.folder):
            logger.warn("A folder named {:s} already exists".format(self.folder))
        else:
            os.makedirs(self.folder, exist_ok=False)
            logger.info("Created folder {:s} for system {:s}".format(self.folder, self.name))
        self.update()

    def set_folder(self, newfolder:str, overwrite = True) -> None:
        if not os.path.isdir(newfolder):
            if os.path.exists(newfolder):
                os.remove(newfolder)
                logger.warn("A file with name {:s} has been removed".format(newfolder))
            os.makedirs(newfolder, exist_ok=False)
        self.folder = newfolder
        self.update(overwrite=True)

    def format_names(self, overwrite = False):
        """ Rename all files in the current system to the system's name """
        for k in self._FILEATTR_FILEPATH_DICT:
            v = self._FILEATTR_FILEPATH_DICT[k]
            dirname = os.path.dirname(v)
            filname = os.path.basename(v)
            extensi = os.path.splitext(filname)[-1]
            newname = os.path.join(dirname, self.name, extensi)
            if newname == v:
                continue
            if os.path.exists(newname):
                if overwrite:
                    os.remove(newname)
                    logger.info("Removed the existing file {}".format(newname))
                    os.rename(v, newname)
                    logger.info("Rename {}.{} from {} to {}".format(self.name, k, v, newname))
                else:
                    logger.info("Skipped the existing file {}".format(newname))
            else:
                os.rename(v, newname)
                logger.info("Rename {}.{} from {} to {}".format(self.name, k, v, newname))

    @property
    def files(self) -> dict:
        res = {}
        for k in self._FILEATTR_FILEPATH_DICT:
            v = self._FILEATTR_FILEPATH_DICT[k]
            if isinstance(v, str) and os.path.isfile(v):
                res[k] = v
        return res
    
    @property
    def reference_name(self) -> str:
        return os.path.join(self.folder, self.name)

# @dataclass 
class Molecule(System): 
    #constants 
    _SUPPORT_INPUT_FORMATS = ['pdb', 'xyz'] 
    # other
    def __init__(
            self, 
            xyzfile:        str, 
            charge:         int, 
            multiplicity:   int,
            name            = "",
            residue_name    = "MOL",
            folder          = WORKING_DIR,
            ) -> None:

        self.name           = process_system_name(name, xyzfile, support_input_format=Molecule._SUPPORT_INPUT_FORMATS)
        self.folder         = os.path.abspath(folder)
        self.charge         = charge
        self.multiplicity   = multiplicity
        self.spinmult       = multiplicity
        self.residue_name   = residue_name
        self.number         = 0
        self.read_coordinate(xyzfile)
        super(Molecule, self).__init__(name = self.name)
        # super(Molecule, self).__post_init__()

    def read_coordinate(self, fname:str):
        ext = os.path.splitext(fname)[-1][1:]
        setattr(self, ext, fname)
        for e in Molecule._SUPPORT_INPUT_FORMATS:
            if e == ext:
                continue
            nname = os.path.splitext(fname)[0] + "." + e
            subprocess.run(f"obabel -i {ext} {fname} -o {e} > {nname}", shell = True)
            setattr(self, e, nname)




if __name__ == "__main__":
    pass