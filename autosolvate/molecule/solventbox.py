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
from .molecule import *



class SolventBox(System):
    # constants
    _SUPPORT_INPUT_FORMATS = ["lib", "off"] 

    def __init__(self, 
            solventbox:     str = "",
            frcmod:         str = "",
            name:           str = "", 
            box_name:       str = "SLVBOX",
            folder:         str = WORKING_DIR,
            amber_solvent:  bool = False,
            ) -> None:
        """
        The data class representing a pre-built solvent box.

        Parameters
        ----------
        solventbox : str
            The path to the solvent box file. The solvent box file can be in either .lib or .off format. Ignored if amber_solvent is True.
        frcmod : str
            The path to the frcmod file. Ignored if amber_solvent is True.
        name : str
            The name of the solvent box. if amber_solvent is True, the name should be one of the following: 'water', 'methanol', 'chloroform', 'nma'.
        box_name : str
            The name of the solvent box in the solvent box file.
        folder : str
            The path to the folder containing the solvent box file and the frcmod file.
        amber_solvent : bool
            A boolean flag indicating whether the solvent box is an AMBER solvent or not. If True, the solventbox and frcmod parameters are ignored.
        """
        self.name           = process_system_name(name, solventbox, support_input_format=SolventBox._SUPPORT_INPUT_FORMATS, check_exist = False if amber_solvent else True)
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
            self.box_name = self.get_box_name(self.solventbox)
            self.solventbox = getattr(self, ext)
        else:
            pass
        super(SolventBox, self).__init__(name = self.name)
        self.logger.name = self.__class__.__name__

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
logging.basicConfig(level = INFO, force = True, handlers=[])