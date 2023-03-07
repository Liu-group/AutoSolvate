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


from logging import DEBUG, INFO, WARN, WARNING, CRITICAL
logging.basicConfig(level = INFO, force = True, handlers=[])

logger              = logging.getLogger(name = "Molecule")
# output_handler      = logging.FileHandler(filename = "log.txt", mode = "a", encoding="utf-8")
output_handler      = logging.StreamHandler()
output_formater     = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
output_handler.setFormatter(output_formater)
if len(logger.handlers) == 0:
    logger.addHandler(output_handler)




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
AMBER_SOLVENT_NAME       =  [AMBER_WATER.name, AMBER_METHANOL.name, AMBER_CHLOROFORM.name, AMBER_NMA.name] 
    
