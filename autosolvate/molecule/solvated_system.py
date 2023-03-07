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
from .molecule import *
from .molecule_complex import *

# Define module logger
from logging import DEBUG, INFO, WARN, CRITICAL
logging.basicConfig(level = INFO, force = True, handlers=[])

logger              = logging.getLogger(name = "SolvatedSystem")
# output_handler      = logging.FileHandler(filename = "log.txt", mode = "a", encoding="utf-8")
output_handler      = logging.StreamHandler()
output_formater     = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
output_handler.setFormatter(output_formater)
if len(logger.handlers) == 0:
    logger.addHandler(output_handler)


class SolvatedSystem(System):
    def __init__(
            self, name:str,
            solute:         System          = None,
            solvent:        System          = None,
            cubesize:       int             = 54,
            closeness:      float           = 2.0,
            solute_number:  int             = 1,
            solvent_number: int             = 1,
            folder:         str             = WORKING_DIR,
            ) -> None:
        self.name           = os.path.basename(os.path.splitext(name)[0])      
        self.folder         = os.path.abspath(folder)
        self.solute         = solute
        self.solvent        = solvent
        self.cubesize       = cubesize
        self.closeness      = closeness
        self.solute_number  = solute_number
        self.solvent_number = solvent_number
        self.solutes        = []
        self.solvents       = []
        self.netcharge      = 0
        self.multiplicity   = 1
        super(SolvatedSystem, self).__init__(name = self.name)

        self.check_solute()
        self.check_solvent()
        self.compute_netcharge()

    def compute_netcharge(self) -> None:
        '''
        compute net charge 
        '''
        for solvent in self.solvents:
            if isinstance(solvent, SolventBox):
                self.netcharge += 0
            elif isinstance(solvent, Molecule):
                self.netcharge += solvent.charge * solvent.number
        for solute in self.solutes:
            if isinstance(solute, Molecule):
                self.netcharge += solute.charge * solute.number
            elif isinstance(solute, MoleculeComplex):
                self.netcharge += solute.netcharge * solute.number
        logger.info("Net charge: {}".format(self.netcharge))

    def compute_multiplicity(self) -> None:
        '''
        compute mutiplicity 
        '''
        not_paired_electrons = 0
        for solute in self.solutes:
            solute:Molecule
            not_paired_electrons += (solute.multiplicity - 1) * solute.number
        self.multiplicity = not_paired_electrons + 1
        logger.info(f"Total multiplicity of the molecule is {self.multiplicity}")

    def check_solute(self) -> None:
        '''
        @TODO 
        check solute 
        '''
        if isinstance(self.solute, Molecule) or isinstance(self.solute, MoleculeComplex):
            logger.info("Use solute: {}".format(self.solute.name))
            if isinstance(self.solute_number, int) and self.solute_number > 0:
                logger.info("Number of solute {} : {}".format(self.solute.name, self.solute_number))
                self.solute.number = self.solute_number
            elif isinstance(self.solute_number, list):
                logger.info("Number of solute {} : {}".format(self.solute.name, self.solute_number))
                self.solute.number = self.solute_number[0]
            else:
                logger.error("Solute number {} for solute {} is invalid".format(self.solute_number, self.solute.name))
                raise ValueError("Solute number is invalid")      
            self.solutes.append(self.solute)
        elif isinstance(self.solute, list):
            if not isinstance(self.solute_number, list) or not len(self.solute_number) == len(self.solute):
                logger.error("Solute number {} is invalid".format(self.solute_number))
                raise ValueError("Solute number is invalid")
            for n, solute in zip(self.solute_number, self.solute):
                logger.info("Use solute: {}".format(solute.name))
                if n <= 0:
                    logger.error("Solute number {} for solute {} is invalid".format(n, solute.name))
                    raise ValueError("Solute number is invalid")
                logger.info("Use solute: {}".format(solute.name))
                logger.info("Number of solute {} : {}".format(solute.name, n))
                solute.number = n
                self.solutes.append(solute)
        else:
            logger.error("Solute {} is invalid".format(self.solute))
            raise ValueError("Solute is invalid")
    
        if len(self.solutes) > 1:
            logger.warn("More than one solute is given. The relative position between solute molecules will be randomly generated")
            logger.warn("If you want to specify the relative position, please generate a molecule using a single xyz/pdb file contains all fragments.")
        
    def check_solvent(self) -> None:
        '''
        @TODO 
        check solvent 
        '''
        if isinstance(self.solvent, list) and len(self.solvent) == 1:
            self.solvent = self.solvent[0]
        if isinstance(self.solvent, SolventBox):
            logger.info("Use solvent: {}".format(self.solvent.name))
            if self.solvent.amber_solvent:
                logger.info("Use AMBER solvent: {}".format(self.solvent.name))
            elif self.solvent.check_exist("off"):
                logger.info("Use existing solvent box: {}".format(self.solvent.off))
            self.solvents.append(self.solvent)
        elif isinstance(self.solvent, Molecule) or isinstance(self.solvent, MoleculeComplex):
            logger.info("Use custom built solvent: {}".format(self.solvent.name))
            self.solvents.append(self.solvent)
            if isinstance(self.solvent_number, int):
                self.solvent.number = self.solvent_number
            elif isinstance(self.solvent_number, list) and len(self.solvent_number) == 1:
                self.solvent.number = self.solvent_number[0]
            else:
                logger.error("Solvent number {} is invalid".format(self.solvent_number))
                raise ValueError("Solvent number is invalid")

        elif isinstance(self.solvent, list) and len(self.solvent) > 1:
            logger.info("More than one solvents used. Will not use pre-built solvent box.")
            for solvent in self.solvent:
                if isinstance(solvent, SolventBox):
                    logger.critical("Cannot use pre-built solvent box with more than one solvents.")
                    raise ValueError("Cannot use pre-built solvent box with more than one solvents.")
            if not isinstance(self.solvent_number, list) or not len(self.solvent_number) == len(self.solvent):
                logger.error("Solvent number {} is invalid".format(self.solvent_number))
                raise ValueError("Solvent number is invalid")
            for n, solvent in zip(self.solvent_number, self.solvent):
                logger.info("Use solvent: {}".format(solvent.name))
                if n <= 0:
                    logger.error("Solvent number {} for solvent {} is invalid".format(n, solvent.name))
                    raise ValueError("Solvent number is invalid")
                logger.info("Number of solvent {} : {}".format(solvent.name, n))
                solvent.number = n
                self.solvents.append(solvent)
        else:
            logger.error("Solvent {} is invalid".format(self.solvent))
            raise ValueError("Solvent is invalid")
        
        if len(self.solutes) > 1 or max([solute.number for solute in self.solutes]) > 1:
            for solvent in self.solvents:
                if isinstance(solvent, SolventBox):
                    logger.critical("Cannot use pre-built solvent box with more than one solute molecule.")
                    raise ValueError("Cannot use pre-built solvent box with more than one solute molecule.")
                
    def has_solute(self, mol:System) -> bool:
        for solute in self.solutes:
            if mol.name == solute.name or id(mol) == id(solute):
                return True
        return False
    
    def has_solvent(self, mol:System) -> bool:
        for solvent in self.solvents:
            if mol.name == solvent.name or id(mol) == id(solvent):
                return True
        return False

    def add_solute(self, mol: System, number: int = 0) -> None: 
        '''
        @TODO 
        check solute validity  
        ''' 
        if self.has_solute(mol):
            logger.warn("Solute {} already exists. Do not add the given solvent.".format(mol.name))
            return 
        if number > 0:
            mol.number = number
        if mol.number <= 0: 
            logger.critical("Solute number {} of solute {} is invalid".format(mol.number, mol.name))
            raise ValueError("Solute number is invalid")

    def add_solvent(self, mol: System, number: int = 0) -> None: 
        '''
        @TODO 
        check solvent validity 
        ''' 
        if isinstance(mol, SolventBox):
            logger.warn("The number argument is ignored because a pre-built solvent box is used.")
            if len(self.solvents) == 0:
                logger.info("Add solventbox: {}".format(mol.name))
                self.solvents.append(mol)
            else:
                logger.critical("Cannot add a solvent box when there is more than one kind of solvent exist.")
                raise ValueError("Cannot add a solvent box when there is more than one kind of solvent exist.")
        elif isinstance(mol, Molecule):
            if self.has_solvent(mol):
                logger.warn("Solvent {} already exists. Do not add the given solvent.".format(mol.name))
                return
            for solvent in self.solvents:
                if isinstance(solvent, SolventBox):
                    logger.critical("Cannot add another solvent because solvent {} is a pre-built solvent box.".format(solvent.name))
                    raise ValueError("Cannot add another solvent because solvent {} is a pre-built solvent box.".format(solvent.name))
        if number > 0 and isinstance(mol, Molecule):
            mol.number = number
        if mol.number <= 0 and isinstance(mol, Molecule):
            logger.critical("Solvent number {} of solvent {} is invalid".format(mol.number, mol.name))
            raise ValueError("Solvent number is invalid")
        
    def set_solute_number(self, solute:System, number: int = 0) -> None:
        if not self.has_solute(solute):
            logger.warn("Solute {} does not exist. Do not set the number.".format(solute.name))
            return
        if number > 0:
            solute.number = number
        else:
            logger.warn("Solute number {} is invalid".format(number))
            return

    def set_solvent_number(self, solvent:System, number: int = 0) -> None:
        if not self.has_solvent(solvent):
            logger.warn("Solvent {} does not exist. Do not set the number.".format(solvent.name))
            return
        if number > 0:
            if isinstance(solvent, SolventBox):
                logger.warn("The using solvent is a pre-built solvent box. The number argument is ignored.")
                return 
            solvent.number = number
        else:
            logger.warn("Solvent number {} is invalid".format(number))
            return

    def set_closeness(self, closeness = 2.0, automate: bool = False) -> None: 
        r'''
        @TODO: 
            1. It only supports one solvent 
        '''     
        if automate: 
            if len(self.solvents) > 1: 
                logger.warn("automatically set closeness only supports one solvent")
                logger.warn("use default closeness: {}".format(self.closeness))
                self.closeness = closeness
            else:
                solvent = self.solvents[0]
                if solvent.name == 'acetonitrile':
                    self.closeness = 1.88 
                elif solvent.name == 'water': 
                    self.closeness = 0.50
                elif solvent.name == 'methanol': 
                    self.closeness = 0.60 
                elif solvent.name == 'nma': 
                    self.closeness = 0.58 
                elif solvent.name == 'chloroform': 
                    self.closeness = 0.58 
                else: 
                    logger.warn("automatically set closeness only supports water, acetonitrile, methanol, nma, chloroform")
                    logger.warn("use default closeness for solvent {}: {}".format(solvent.name, self.closeness))
        else: 
            logger.info("set closeness manually")
            logger.info("use closeness: {}".format(closeness))
            self.closeness = closeness
    
    # def update(self) -> None:
    # 这个方法是干嘛的？我觉得没用，就删了。（这句话又是 copilot 乱加的。）
    # The update method is already defined in the parent class.