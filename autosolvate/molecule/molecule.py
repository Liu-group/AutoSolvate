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
            self, xyzfile:str, charge:int, multiplicity:int,
            name = "",
            residue_name = "MOL",
            folder = WORKING_DIR,
            ) -> None:
        if os.path.isfile(xyzfile) and os.path.splitext(xyzfile)[-1][1:] in Molecule._SUPPORT_INPUT_FORMATS:
            self.read_coordinate(name)
        else:
            raise ValueError("The input file {:s} does not exist".format(xyzfile))
        if name == "":
            name = os.path.basename(os.path.splitext(name)[0])  
        self.name           = name      
        self.folder         = os.path.abspath(folder)
        self.charge         = charge
        self.multiplicity   = multiplicity
        self.spinmult       = multiplicity
        self.residue_name   = residue_name
        self.number         = 0
        self.is_solventbox  = False
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


class SolventBox(System):
    # constants
    _SUPPORT_INPUT_FORMATS = ["lib", "off"] 

    def __init__(
            self, name:str, solventbox:str,
            frcmod:         str = "",
            box_name:       str = "SLVBOX",
            folder:         str = WORKING_DIR,
            amber_solvent:  bool = False,
            ) -> None:
        self.name           = os.path.basename(os.path.splitext(name)[0])      
        self.folder         = os.path.abspath(folder)
        self.frcmod         = os.path.abspath(frcmod)
        self.box_name       = box_name
        self.amber_solvent  = amber_solvent
        self.solventbox     = os.path.abspath(solventbox)
        

        if not amber_solvent:
            if not os.path.exists(solventbox):
                logger.error("The pre-built solvent box file {} does not exist".format(solventbox))
                raise FileNotFoundError("The solvent box file does not exist")
            ext = os.path.splitext(solventbox)[-1][1:]
            if ext not in SolventBox._SUPPORT_INPUT_FORMATS:
                logger.error("The input file format of the solvent box is not supported")
                raise ValueError("The input file format of the solvent box is not supported")
            setattr(self, ext, solventbox)
            # self.reset_box_name(self.solventbox)
            self.box_name = self.get_box_name(self.solventbox)
            self.solventbox = getattr(self, ext)
        else:
            pass
        super(SolventBox, self).__init__(name = self.name)
        # super(SolventBox, self).__post_init__()

    def get_box_name(self, fname:str):
        with open(fname, "r") as f:
            lines = f.readlines()
        rline, nline = 0, 0
        for line in lines:
            if line.startswith("!"):
                rline += 1
            else:
                nline += 1
                rline += 1
            if nline == 1:
                logger.info("the name of the solvent box is found at line {}".format(rline))
                previousname = line.strip().replace('"', "")
                logger.info("The previous name of the solvent box is {}".format(previousname))
                break
        return previousname

    def reset_box_name(self, fname:str):
        with open(fname, "r") as f:
            lines = f.readlines()
        rline, nline = 0, 0
        for line in lines:
            if line.startswith("!"):
                rline += 1
            else:
                nline += 1
                rline += 1
            if nline == 1:
                logger.info("the name of the solvent box is found at line {}".format(rline))
                previousname = line.strip().replace('"', "")
                logger.info("The previous name of the solvent box is {}".format(previousname))
                break
        logger.info("The new name of the solvent box is {}".format(self.box_name))
        with open(fname, "w") as f:
            for i in range(0, rline - 1):
                f.write(lines[i])
            f.write('"{}"\n'.format(self.box_name))
            for i in range(rline, len(lines)):
                print(lines[i].startswith("!entry.{}".format(previousname)), "!entry.{}".format(previousname), lines[i])
                if lines[i].startswith("!entry.{}".format(previousname)):
                    newline = lines[i].replace("!entry.{}".format(previousname), "!entry.{}".format(self.box_name))
                else:
                    newline = lines[i]
                f.write(newline)
        logger.info("The solvent box file {} has been updated".format(fname))



AMBER_WATER               = SolventBox( name='water',
                                        solventbox='solvents.lib',
                                        box_name='TIP3PBOX',
                                        amber_solvent = True
                                    )

AMBER_METHANOL            = SolventBox( name='methanol',
                                        solventbox='solvents.lib',
                                        frcmod="frcmod.meoh",
                                        box_name='MEOHBOX',
                                        amber_solvent = True
                                    )                           

AMBER_CHLOROFORM          = SolventBox(name='chloroform',
                                        solventbox='solvents.lib',
                                        frcmod="frcmod.chcl3",
                                        box_name='CHCL3BOX',
                                        amber_solvent = True
                                    )
                            
AMBER_NMA                 = SolventBox(name='nma',
                                        solventbox='solvents.lib',
                                        frcmod="frcmod.nma",
                                        box_name='NMABOX',
                                        amber_solvent = True
                                    )

AMBER_SOLVENT_LIST       =  [AMBER_WATER, AMBER_METHANOL, AMBER_CHLOROFORM, AMBER_NMA] 
AMBER_SOLVENT_DICT       =  {
                             AMBER_WATER.name       :AMBER_WATER, 
                             AMBER_METHANOL.name    :AMBER_METHANOL, 
                             AMBER_CHLOROFORM.name  :AMBER_CHLOROFORM, 
                             AMBER_NMA.name         :AMBER_NMA,
                            }
    

if __name__ == "__main__":
    pass