# This module is here to process the json input of multicomponent.py

import os 
import re
import copy
import json
import shutil 
import logging

from typing import Iterable, List, Dict, Any, Tuple

from .tools import *
from ..Common import *


def convert_cmd_to_dict(argument_dict:Dict[str, Any]) -> Dict[str, Any]:
    special_arguments = ["solute", "solvent"]
    newdata = {}
    for key, value in argument_dict.items():
        if key not in special_arguments:
            newdata[key] = value
    data = argument_dict
    if "main" in data and data["main"]:
        newdata["solute"] = dict()
        newdata['solute']['xyzfile'] = data["main"]
        if "charge" in data:
            newdata["solute"]["charge"] = data["charge"]
        if "spinmult" in data:
            newdata["solute"]["spinmult"] = data["spinmult"]
    if "solvent" in data and data["solvent"]:
        newdata["solvent"] = dict()
        if "solventoff" in data:
            newdata["solvent"]["off"] = data["solventoff"]
        if "solventfrcmod" in data:
            newdata["solvent"]["frcmod"] = data["solventfrcmod"]
        if "solvent" in data:
            newdata["solvent"]["name"] = data["solvent"]
    return newdata


def determine_mw_from_xyz(xyzfile:str) -> float:
    """
    Determine the molecular weight of the molecule from the xyz file

    Parameters
    ----------
    xyzfile : str
        xyz file name

    Returns
    -------
    mw : float
        molecular weight in g/mol
    """
    with open(xyzfile, "r") as f:
        lines = f.readlines()
    mw = 0
    for line in lines[2:]:
        parts = line.split()
        if len(parts) == 4:
            if parts[0] not in ATOMIC_MW:
                raise ValueError(f"Cannot determine the molecular weight of the atom {parts[0]}")
            mw += ATOMIC_MW[parts[0]]
    return mw


def add_missing_xyzfile_keyword(data:dict, support_input_format:Iterable[str] = ("xyz", "pdb", "mol2", "prep", "off", "lib")) -> dict:
    if "xyzfile" in data:
        return data
    for key, value in data.items():
        if key in support_input_format and not "xyzfile" in data and isinstance(value, str) and os.path.isfile(value):
            data["xyzfile"] = value
            break
    return data

def overwrite_keyword(data:dict, arguments:List[Tuple[str, Any]]) -> dict:
    """
    Overwrite the keyword in the dictionary with the arguments
    """
    newdata = copy.deepcopy(data)
    ignoredargs = []
    arguments
    for key, value in arguments:
        if isinstance(value, str) and value == "":
            ignoredargs.append(key)
    for key in ignoredargs:
        arguments.remove((key, ""))
    for key, value in arguments:
        newdata[key] = value
    if "main" in data and data["main"]:
        newdata["solute"] = dict()
        newdata['solute']['xyzfile'] = data["main"]
        if "charge" in data:
            newdata["solute"]["charge"] = data["charge"]
        if "spinmult" in data:
            newdata["solute"]["spinmult"] = data["spinmult"]
    if "solvent" in data and data["solvent"]:
        newdata["solvent"] = dict()
        if "solventoff" in data:
            newdata["solvent"]["off"] = data["solventoff"]
        if "solventfrcmod" in data:
            newdata["solvent"]["frcmod"] = data["solventfrcmod"]
        if "solvent" in data:
            newdata["solvent"]["name"] = data["solvent"]

    return newdata

def check_inputs(data:dict) -> dict:
    """
    Check if the input is valid
    """
    if "solute" not in data and "solutes" not in data:
        raise ValueError("No solute is provided")
    if "solvent" not in data and "solvents" not in data:
        raise ValueError("No solvent is provided")
    if "solute" in data and "solutes" in data:
        raise ValueError("You cannot provide both solute and solutes, only provide the solute argument if you want to use single solute")
    if "solvent" in data and "solvents" in data:
        raise ValueError("You cannot provide both solvent and solvents, only provide the solvent argument if you want to use single solvent")
    if "solutes" in data and (not isinstance(data["solutes"], list) or len(data["solutes"]) == 0):
        raise ValueError("solutes must be a list")
    if "solvents" in data and (not isinstance(data["solvents"], list) or len(data["solvents"]) == 0):
        raise ValueError("solvents must be a list")
    


class InputParser(object):

    def __init__(self):
        

        self.data = {}
        self.logger = logging.getLogger(name = self.__class__.__name__)
        self.output_handler             = logging.FileHandler(filename = "autosolvate.log", mode = "a", encoding="utf-8")
        self.output_formater            = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
        self.output_handler.setFormatter(self.output_formater)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(self.output_handler)

        self.warn_user = False

    def correct_keyword(self, data: dict) -> dict:
        """
        Recursively correct the keys in the dictionary and its subdictionaries to the acceptable keys
        """
        def correct_dict(d: dict) -> dict:
            newdata = {}
            for key, value in d.items():
                if isinstance(value, dict):
                    value = correct_dict(value)
                if key in keyword_dict:
                    newdata[keyword_dict[key]] = value
                    self.logger.debug(f"Correcting keyword {key} to {keyword_dict[key]}")
                else:
                    newdata[key] = value
            return newdata

        keyword_dict = {
            "chargemethod": "charge_method",
            "cubesize": "cube_size",
            "fragmentcharge": "fragment_charge",
            "fragmentspin": "fragment_spinmult",
            "fragmentspinmultiplicity": "fragment_spinmult",
            "fragmentmultiplicity": "fragment_spinmult",
            "spinmultiplicity": "spinmult",
            "multiplicity": "spinmult",
            "spin": "spinmult",
            "solventoff": "solvent_off",
            "solventfrcmod": "solvent_frcmod",
            "rundir": "folder",
            "output": "prefix",
            "basis_set": "basisset",
            "qmdir": "QMexe",
            "qmexe": "software",
        }

        return correct_dict(data)


    def read_json(self, file:str):
        with open(file, "r") as f:
            data = json.load(f)
        self.data = data

    def read_dict(self, data:dict):
        self.data = data

    def determine_qm_exe_location(self, qm_software_name:str):
        """
        Determine the location of the QM software executable

        Parameters
        ----------
        qm_software_name : str
            name of the QM software, e.g. gau,g09,g16,gms,orca
        """
        if qm_software_name.lower() in ["gaussian", "gau" ]:
            qm_software_name = "gau"
        elif qm_software_name.lower() in ["gamess-us", "gamess", "rungms"]:
            qm_software_name = "rungms"
        else:
            qm_software_name = qm_software_name.lower()
        
        if qm_software_name in ["gamess-us", "rungms", "gamess"]:
            qm_software_name_ = "GAMESS"
        else:
            qm_software_name_ = qm_software_name
        # use which command to determine the location of the executable
        if not shutil.which(qm_software_name):
            # try to load the software with the "module load" command
            if shutil.which("module"):
                os.system(f"module load {qm_software_name_}")
        if shutil.which(qm_software_name):
            self.logger.info(f"Found the executable of {qm_software_name} at {shutil.which(qm_software_name)}")
            return shutil.which(qm_software_name)
        else:
            self.logger.error(f"Cannot find the executable of {qm_software_name}. Please provide the path of the executable.")
            raise FileNotFoundError(f"Cannot find the executable of {qm_software_name}")
        
    def add_qm_kwargs(self):
        default_qm_kwargs = {
            "method": "b3lyp",
            "basisset": "DEF2-TZVP",
            "maxcore": 1024,
            "nprocs": 1,
            "opt": True,
        }
        if 'software' not in self.data:
            self.warn_user = True
            self.data['software'] = "orca"
        if 'QMexe' not in self.data or not os.path.exists(self.data['QMexe']):
            path = self.determine_qm_exe_location(self.data['software'])
            self.data['QMexe'] = path
        not_all_parameters = False
        for key, value in default_qm_kwargs.items():
            if key not in self.data:
                not_all_parameters = True
                self.data[key] = value
        if not_all_parameters:
            self.warn_user = True
        self.data["qm_kwargs"] = {
            "basisset": self.data["basisset"],
            "method": self.data["method"],
            "software": self.data["software"],
            "QMexe": self.data["QMexe"],
            "maxcore": self.data["maxcore"],
            "nprocs": self.data["nprocs"],
            "opt": self.data["opt"],
        }

    def add_mcpb_kwargs(self):
        default_mcpb_kwargs = {
            "cutoff": 2.8,
            "fakecharge": False,
            "mode": "A",
        }
        # check if amberhome exists 
        if 'amberhome' not in self.data:
            self.data['amberhome'] = os.getenv("AMBERHOME") + "/bin/"
        amberhomepath = os.path.expandvars(self.data['amberhome'])
        if not os.path.exists(amberhomepath):
            self.logger.error(f"Cannot find the amberhome path {amberhomepath}. Make sure AMBER is installed and the path is correct.")
            raise FileNotFoundError(f"Cannot find the amberhome path {amberhomepath}")
        not_all_parameters = False
        for key, value in default_mcpb_kwargs.items():
            if key not in self.data:
                not_all_parameters = True
                self.data[key] = value
        if not_all_parameters:
            self.warn_user = True
        self.data["mcpb_kwargs"] = {
            "amberhome": amberhomepath,
            "cutoff": self.data["cutoff"],
            "fakecharge": self.data["fakecharge"],
            "mode": self.data["mode"],
        }

    def add_tmc_kwargs(self, solutedata:dict):
        if "metal_charge" not in solutedata:
            raise ValueError("The 'metal_charge' paremeter is required to set the valance of the metal atom.")
        if "spinmult" not in solutedata:
            raise ValueError("The 'spinmult' or 'spin' parameter is required to set the spin multiplicity of the transition metal complex.")
        if "total_charge" not in solutedata:
            solutedata["total_charge"] = "default"
            self.logger.info("The 'total_charge' parameter is not provided. Will be determined automatically.")
        if "chargefile" not in solutedata:
            solutedata["chargefile"] = ""    

    def assign_solvent_numbers(self):
        # in-place operation
        all_number_exist = True
        all_number_missing = True
        all_vratio_exist = True
        all_wratio_exist = True
        for solvent in self.data['solvents']:
            if "number" not in solvent:
                all_number_exist = False
            if "number" in solvent:
                all_number_missing = False
            if "volume_ratio" not in solvent:
                all_vratio_exist = False
            if "weight_ratio" not in solvent:
                all_wratio_exist = False
            if "name" not in solvent:
                if "xyzfile" in solvent:
                    solvent["name"] = os.path.basename(solvent["xyzfile"]).split(".")[0]
                else:
                    raise ValueError("Both 'name' and 'xyzfile' is missing for a solvent.")

        if not all_number_exist and not all_number_missing:
            raise ValueError("Some solvents have the 'number' parameter while others do not. Please provide the 'number' parameter for all solvents.")
        if all_number_exist:
            return
        # all solvent numbers are missing. compute from the weight ratio
        if "cube_size" not in self.data:
            raise ValueError("The 'cubesize' parameter is required to compute the number of solvent molecules.")
        if not isinstance(self.data["cube_size"], Iterable):
            volume = float(self.data["cube_size"]) ** 3
        else:
            a, b, c = self.data["cube_size"]
            volume = float(a) * float(b) * float(c)
        volume_in_m3 = volume * 1e-30
        if len(self.data['solvents']) == 1:
            solventname = self.data['solvents'][0]['name']
            if solventname in AMBER_SOLVENT_NAMES:
                return  # Use amber prebuilt solvent box.
            elif solventname in SOLVENT_DENSITY:
                number  = calculate_solvent_number(solventname, volume_in_m3)
                self.data['solvents'][0]['number'] = number
                return 
            else:
                raise ValueError(f"Cannot determine the number of the solvent. Please provide the 'number' parameter.")
        # if there are multiple solvents, the number of each solvent must be provided
        if not all_vratio_exist and not all_wratio_exist:
            raise ValueError("Please provide the 'volume_ratio' or 'weight_ratio' parameter for all solvents if the 'number' parameter is not provided.")
        for solvent in self.data['solvents']:
            # determine density
            if solvent["name"] in SOLVENT_DENSITY and "density" not in solvent:
                solvent["density"] = SOLVENT_DENSITY[solvent["name"]]
            elif "density" not in solvent:
                raise ValueError(f"The density of the solvent {solvent['name']} is not provided.")
            # determine the molecular weight. 
            if solvent["name"] in SOLVENT_MW and "molecular_weight" not in solvent:
                solvent["molecular_weight"] = SOLVENT_MW[solvent["name"]]
            elif "molecular_weight" in solvent:
                pass    # sometimes a solvent may contain isotopes, so the user may provide the molecular weight
            elif "molecular_weight" not in solvent and "xyzfile" in solvent:
                solvent["molecular_weight"] = determine_mw_from_xyz(solvent["xyzfile"])
            else:
                raise ValueError(f"The molecular weight of the solvent {solvent['name']} is not provided.")
        # determine the number of solvent molecules
        if all_vratio_exist:
            numbers = calculate_solvent_numbers_from_volume_portions(
                    [s["molecular_weight"] for s in self.data["solvents"]],
                    [s["density"] for s in self.data["solvents"]],
                    [s["volume_ratio"] for s in self.data["solvents"]],
                    volume_in_m3
                )
        elif all_wratio_exist:
            numbers = calculate_solvent_numbers_from_weight_portions(
                    [s["molecular_weight"] for s in self.data["solvents"]],
                    [s["density"] for s in self.data["solvents"]],
                    [s["weight_ratio"] for s in self.data["solvents"]],
                    volume_in_m3
                )
        for i, solvent in enumerate(self.data['solvents']):
            solvent['number'] = numbers[i]

    def parse(self):
        # step 1: correct the keyword
        self.data = self.correct_keyword(self.data)

        # step 2: change 'solute': {...} to 'solutes': [{...},]; change 'solvent': {...} to 'solvents': [{...},]
        check_inputs(self.data)
        if ('solute' not in self.data or len(self.data['solute']) == 0) and ('solutes' not in self.data or len(self.data['solutes']) == 0):
            raise ValueError("No solute is provided.")
        if ('solute' in self.data and len(self.data['solute']) > 0) and ('solutes' in self.data and len(self.data['solutes']) > 0):
            raise ValueError("Both 'solute' and 'solutes' are provided. Please provide only one.")
        if ('solute' in self.data and len(self.data['solute']) > 0) and ('solutes' not in self.data or len(self.data['solutes']) == 0):
            self.data['solutes'] = [self.data['solute']]
            self.data.pop('solute')

        if ('solvent' not in self.data or len(self.data['solvent']) == 0) and ('solvents' not in self.data or len(self.data['solvents']) == 0):
            raise ValueError("No solvent is provided.")
        if ('solvent' in self.data and len(self.data['solvent']) > 0) and ('solvents' in self.data and len(self.data['solvents']) > 0):
            raise ValueError("Both 'solvent' and 'solvents' are provided. Please provide only one.")
        if ('solvent' in self.data and len(self.data['solvent']) > 0) and ('solvents' not in self.data or len(self.data['solvents']) == 0):
            self.data['solvents'] = [self.data['solvent']]
            self.data.pop('solvent')

        # step 3: add xyzfile keyword if not provided
        self.data['solvents'] = [add_missing_xyzfile_keyword(solvent) for solvent in self.data['solvents']]
        self.data['solutes'] = [add_missing_xyzfile_keyword(solute) for solute in self.data['solutes']]

        # step 4: check solute type. add "__TYPE__" keyword to indicate the type of the solute
        has_tmc = False
        for solute in self.data['solutes']:
            if check_transition_metal_complex(solute['xyzfile']):
                solute["__TYPE__"] = "transition_metal_complex"
                has_tmc = True
            elif check_multicomponent(solute['xyzfile']):
                solute["__TYPE__"] = "complex"
            else:
                solute["__TYPE__"] = "molecule"
        # for solvent in self.data['solvents']:
        #     if check_transition_metal_complex(solvent['xyzfile']) or check_multicomponent(solvent['xyzfile']):
        #         raise ValueError("Transition metal complex or molecule complex as solvent is not supported.")
            
        # step 5: add the missing keywords for transiton metal complex
        if has_tmc:
            self.logger.info("Transition metal complex is detected. Will use AutoMCPB to generate the force field parameters.")
        for solute in self.data['solutes']:
            if "number" not in solute:
                solute["number"] = 1
            if solute["__TYPE__"] != "transition_metal_complex":
                continue    # the other type of solvents can be handled by the default parameters
            self.add_tmc_kwargs(solute) # this is an in-place operation on the dictionary, no need to return the value

        # step 6: if there is a transition metal complex, check the parameters required for automcpb.
        if has_tmc:
            self.add_qm_kwargs()
            self.add_mcpb_kwargs()

        # step 7: assign the number of solvent molecules
        self.assign_solvent_numbers()

        # step 8: notify the user if some parameters are missing
        if self.warn_user:
            self.logger.warning("Some parameters are set to default ones that may not work for your system.")
            self.logger.warning("Please check the 'autosolvate_input_full.json' file to see the full list of parameters.")