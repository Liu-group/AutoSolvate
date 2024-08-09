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
import logging

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
        self.logger                     = logging.getLogger(name = "System")
        self.output_handler             = logging.FileHandler(filename = "autosolvate.log", mode = "a", encoding="utf-8")
        # self.output_handler             = logging.StreamHandler()
        self.output_formater            = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
        self.output_handler.setFormatter(self.output_formater)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(self.output_handler)

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
        if isinstance(__value, logging.Logger):
            pass
        elif isinstance(__value, str):
            if os.path.isfile(__value) and len(__name) <= 6:
                __value = os.path.abspath(__value)
                self._FILEATTR_FILEPATH_DICT[__name] = __value
                self.logger.debug("set the '{:s}' file of system '{:s}' to '{:s}'".format(__name, self.name, __value))
                if not self.check_extension(__name, __value):
                    self.logger.warning("The '{:s}' file of system '{:s}' does not have the correct extension".format(__name, self.name, __value))
            else:
                if __name in self._FILEATTR_FILEPATH_DICT:
                    self._FILEATTR_FILEPATH_DICT[__name] = __value
                    if not (hasattr(self, "amber_solvent") and self.amber_solvent):
                        self.logger.warning("The '{:s}' file of system '{:s}' is set to a non-existent path '{:s}'".format(__name, self.name, __value))
                else:
                    self.logger.debug("set the '{:s}' attribute of system '{:s}' as {:s}".format(__name, self.name, __value))
        else:
            if __name in self._FILEATTR_FILEPATH_DICT:
                self.logger.warning("set the '{:s}' file of system '{:s}' as a {:s}".format(__name, self.name, type(__value).__name__))
            else:
                self.logger.debug("set the '{:s}' attribute of system '{:s}' as a {:s}".format(__name, self.name, type(__value).__name__))
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
                self.logger.debug("The '{:s}' file of system '{:s}' is not defined".format(ext, self.name))
            elif not isinstance(self.__getattribute__(ext), str):
                flag *= False
                self.logger.debug("The '{:s}' attribute of system '{:s}' does not exist".format(ext, self.name))
            elif not os.path.isfile(self.__getattribute__(ext)):
                flag *= False
                self.logger.debug("The '{:s}' file of system '{:s}' does not exist".format(ext, self.name))
            else:
                self.logger.debug("The '{:s}' file of system '{:s}' is '{:s}'".format(ext, self.name, self.__getattribute__(ext)))
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
                    self.logger.debug("The '{:s}' file of system '{:s}' already exists".format(ext, self.name))
                elif overwrite:
                    shutil.copy(srcpath, self.folder)
                    self.logger.debug("Copied '{:s}' file of system '{:s}' to '{:s}'".format(ext, srcpath, newpath))
                    self.logger.info("File {:s} was rewritten".format(newpath))
                    self.__setattr__(ext, newpath)
                else:
                    self.logger.debug("The '{:s}' file of system '{:s}' already exists at {:s} ".format(ext, self.__getattribute__(ext), newpath))
                    self.__setattr__(ext, newpath)
            else:
                shutil.copy(srcpath, self.folder)
                self.logger.debug("Copied '{:s}' file of system '{:s}' to '{:s}'".format(ext, srcpath, newpath))
                self.__setattr__(ext, newpath)

    def copy2folder(self, overwrite = True):
        self.update(overwrite=overwrite)

    def _init_folder(self) -> None:
        if self.folder == WORKING_DIR: 
            self.logger.info("Use the current working directory {:s}".format(self.folder))
        elif os.path.exists(self.folder):
            self.logger.debug("A folder named {:s} already exists".format(self.folder))
        else:
            os.makedirs(self.folder, exist_ok=False)
            self.logger.info("Created folder {:s} for system {:s}".format(self.folder, self.name))
        self.update()

    def set_folder(self, newfolder:str, overwrite = True) -> None:
        if not os.path.isdir(newfolder):
            if os.path.exists(newfolder):
                os.remove(newfolder)
                self.logger.warning("A file with name {:s} has been removed".format(newfolder))
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
                    self.logger.info("Removed the existing file {}".format(newname))
                    os.rename(v, newname)
                    self.logger.info("Rename {}.{} from {} to {}".format(self.name, k, v, newname))
                else:
                    self.logger.info("Skipped the existing file {}".format(newname))
            else:
                os.rename(v, newname)
                self.logger.info("Rename {}.{} from {} to {}".format(self.name, k, v, newname))

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
    _SUPPORT_INPUT_FORMATS = ['pdb', 'xyz', 'mol2', 'prep']
    # other
    def __init__(
            self, 
            xyzfile:        str, 
            charge:         int, 
            multiplicity:   int,
            name            = "",
            residue_name    = "",
            folder          = WORKING_DIR,
            amber_solvent   = False,
            ) -> None:
        """
        The data class for holding all files of a molecule.

        Parameters
        ----------
        xyzfile : str
            The file path of the structure file of the molecule, can be in xyz, pdb, mol2 or prep format.
        charge : int
            The charge of the molecule.
        multiplicity : int
            The multiplicity of the molecule.
        name : str, optional
            The name of the molecule, by default "" will be the basename of the xyzfile.
        residue_name : str, optional
            The residue name of the molecule, ignored if the xyzfile is provided as mol2 or prep format, by default will be assigned according to the filename.
        folder : str, optional
            The folder to store all files of the molecule, by default is the current working directory.
        amber_solvent : bool, optional
            A flag to indicate whether the molecule will use amber predefined parameters, by default False.
        """

        self.folder         = os.path.abspath(folder)
        self.charge         = charge
        self.multiplicity   = multiplicity
        self.spinmult       = multiplicity
        self.residue_name   = residue_name
        self.number         = 0
        self.amber_solvent = amber_solvent 
        self.name           = name
        super(Molecule, self).__init__(name = self.name)
        self.logger.name = self.__class__.__name__

        if not self.amber_solvent:
            if not self.name:
                self.name = process_system_name(name, xyzfile, support_input_format=Molecule._SUPPORT_INPUT_FORMATS, check_exist=False)  
            self.process_input_xyz(xyzfile)
            self.generate_pdb()

    def process_input_xyz(self, xyzfile:str) -> None:
        for e in Molecule._SUPPORT_INPUT_FORMATS:
            ext = os.path.splitext(xyzfile)[-1]
            if e in ext:
                self.__setattr__(e, xyzfile)
                break
        else:
            raise ValueError(f"Unsupported input format for the molecule: {xyzfile}")

    def generate_pdb(self):
        """
        Convert the input structure file to pdb file if the provided coordinate file is not in pdb format.
        If the pdb file for the same molecule already exists, it will be ignored.
        This behavior is to keep the residue name of the molecule consistent.
        """
        if self.name == "water" and self.amber_solvent:
            assign_water_pdb(self)
            self.logger.info(f"AMBER predefined water is used.")
            self.logger.info(f"Write the reference pdb file for {self.name} to {self.pdb}")
            self.logger.info(f"Water model used: TIP3P")
        elif self.amber_solvent:
            prep2pdb(self)
            self.logger.info(f"Predefined solvent {self.name} is used.")
            self.logger.info(f"Converted the predefined prep file {self.prep} to pdb file {self.pdb}")
        elif self.check_exist("prep"):
            newpath = self.reference_name + ".pdb"
            if self.check_exist("pdb"):
                self.logger.warning(f"The existing pdb file {self.pdb} will be ignored as amber prep file is provided.")
                shutil.move(self.pdb, self.reference_name + "-bak.pdb")
            self.logger.info(f"Converted the prep file {self.prep} to pdb file {self.pdb}")
            prep2pdb(self)
        elif self.check_exist("lib") or self.check_exist("off"):
            newpath = self.reference_name + ".pdb"
            if self.check_exist("pdb"):
                self.logger.warning(f"The existing pdb file {self.pdb} will be ignored as amber lib/off file is provided.")
                shutil.move(self.pdb, self.reference_name + "-bak.pdb")
            self.logger.info(f"Converted the lib file {self.lib} to pdb file {self.pdb}")
            lib2pdb(self)
        elif self.check_exist("mol2"):
            newpath = self.reference_name + ".pdb"
            if self.check_exist("pdb"):
                self.logger.warning(f"The existing pdb file {self.pdb} will be ignored as mol2 file is provided.")
                shutil.move(self.pdb, self.reference_name + "-bak.pdb")
            self.logger.info(f"Converted the mol2 file {self.mol2} to pdb file {self.pdb}")
            subprocess.run(f"obabel -i mol2 {self.mol2} -o pdb -O {newpath} ---errorlevel 0", shell = True)
            self.pdb = newpath
        elif self.check_exist("pdb"):
            pass
        elif self.check_exist("xyz"):
            newpath = self.reference_name + ".pdb"
            subprocess.run(f"obabel -i xyz {self.xyz} -o pdb -O {newpath} ---errorlevel 0", shell = True)
            self.pdb = newpath

    def get_residue_name(self):
        """
        Get the residue name from the provided pdb, mol2, or prep file.
        Note that the residue name of the molecule will be automatically assigned if the residue name is not found.
        Note if the mol2, or prep file is provided, the 'residue_name' attribute will be updated according to these files instead of the argument 'residue_name'.
        """
        new_residue_name = ""
        if self.check_exist("prep"):
            new_residue_name = extract_residue_name_from_prep(self.prep)
        elif self.check_exist("lib"):
            new_residue_name = extract_residue_name_from_lib(self.lib)
        elif self.check_exist("off"):
            new_residue_name = extract_residue_name_from_lib(self.off)
        elif self.check_exist("mol2"):
            newpath = self.reference_name + "-frommol2.pdb"
            subprocess.run(f"obabel -i mol2 {self.mol2} -o pdb -O {newpath} ---errorlevel 0", shell = True)
            new_residue_name = get_residue_name_from_pdb(newpath)
        elif self.residue_name:
            new_residue_name = self.residue_name
        else:
            self.logger.warning(f"Cannot find the residue name of molecule {self.name}, trying to automatically assign one.")
            new_residue_name = try_ones_best_to_get_residue_name(self, self.pdb, self.name)
        if not new_residue_name:
            new_residue_name = "MOL"
        if self.residue_name != new_residue_name:
            self.logger.warning(f"Residue name of {self.name} changed: {self.residue_name} -> {new_residue_name}")
            self.logger.warning(f"This behavior is to keep the residue name of the molecule consistent.")
            self.residue_name = new_residue_name
        self.logger.info(f"Residue name of {self.name} is set to {self.residue_name}")

    def update(self):
        if not self.amber_solvent:
            super(Molecule, self).update()



# if one want to use other water forcefield, one may need to change the residue name of water to other, such as PL3 for POL3 model. 
AMBER_WATER             = Molecule(name = "water",
                                   xyzfile = None,
                                   charge = 0,
                                   multiplicity = 1,
                                   residue_name = "WAT",
                                   folder = WORKING_DIR,
                                   amber_solvent = True, 
                                  )

AMBER_METHANOL          = Molecule(name = "methanol",
                                   xyzfile = None, 
                                   charge = 0,
                                   multiplicity = 1,
                                   residue_name = "MOH",
                                   folder = WORKING_DIR,
                                   amber_solvent = True, 
                                  )
AMBER_METHANOL.frcmod = "frcmod.meoh"
AMBER_METHANOL.prep = "meoh.in"

AMBER_CHLOROFORM        = Molecule(name = "chloroform",
                                   xyzfile = None,
                                   charge = 0,
                                   multiplicity = 1,
                                   residue_name = "CL3",
                                   folder = WORKING_DIR,
                                   amber_solvent = True, 
                                  )
AMBER_CHLOROFORM.frcmod = "frcmod.chcl3"
AMBER_CHLOROFORM.prep = "chcl3.in"

AMBER_NMA              = Molecule(name = "nma",
                                   xyzfile = None,
                                   charge = 0,
                                   multiplicity = 1,
                                   residue_name = "NMA",
                                   folder = WORKING_DIR,
                                   amber_solvent = True, 
                                  )
AMBER_NMA.frcmod = "frcmod.nma"
AMBER_NMA.prep = "nma.in"

AMBER_SOLVENT_LIST        = [AMBER_WATER, AMBER_METHANOL, AMBER_CHLOROFORM, AMBER_NMA]  
AMBER_SOLVENT_DICT        = { 
                            AMBER_WATER.name:       AMBER_WATER,
                            AMBER_METHANOL.name:    AMBER_METHANOL,
                            AMBER_CHLOROFORM.name:  AMBER_CHLOROFORM,
                            AMBER_NMA.name:         AMBER_NMA,
                            }

if __name__ == "__main__":
    pass