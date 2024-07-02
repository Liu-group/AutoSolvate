import os 
import re
import shutil 


def parse_multisolvent_input(filepath:str):
    """
    Parse the input file for defining the multisolvent system
    
    Parameters
    ----------
    filepath : str
        The path to the input file.
    
    Returns
    -------
    dict
        The dictionary containing the parsed input.
    """
    solvent_input = {}
    with open(filepath, "r") as f:
        for line in f:
            args = line.strip().split()
            solvent_dict = {}
            solvent_dict["xyzpath"] = args[0]
            solvent_dict["number"] = int(args[1])
            solvent_dict[""]
            solvent_dict["charge"] = int(args[2])
            solvent_dict["spinmult"] = int(args[3])