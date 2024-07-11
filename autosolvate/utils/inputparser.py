# This module is here to process the json input of multicomponent.py

import os 
import re
import copy
import json
import shutil 

from typing import Iterable, List, Dict, Any, Tuple


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


def correct_keyword(data:dict) -> dict:
    """
    Correct the key in the dictionary to the acceptable keys
    """
    keyword_dict = {
        "chargemethod": "charge_method",
        "cubesize": "cube_size",
        "fragmentcharge": "fragment_charge",
        "fragmentspinmultiplicity": "fragment_spinmult", 
        "fragmentmultiplicity": "fragment_spinmult",
        "spinmultiplicity": "spinmult",
        "multiplicity": "spinmult",
        "output": "prefix", 
    }
    newdata = {}
    for key, value in data.items():
        if key in keyword_dict:
            newdata[keyword_dict[key]] = value
        else:
            newdata[key] = value
    return newdata

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