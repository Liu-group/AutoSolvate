#@TODO: 
#       is a better way to implement this? : central control of the whole project
#       1. DRY_RUN: is not working 
import os 


#global variable for the whole project 
DRY_RUN         = False
USE_SRUN        = False

WORKING_DIR     = os.getcwd() + '/'

AMBER_SOLVENT_NAMES = ["water", "methanol", "chloroform", "nma"]
AMINO_ACID_RESIDUES = set(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"])

# external programs
# ambertools
ANTECHAMBER     = "antechamber"
PACKMOL         = "packmol"
TLEAP           = "tleap"
PARMCHK         = "parmchk2"
MDGX            = "mdgx"

# QM programs
GAUSSIAN        = "g16"
TERACHEM        = "terachem"

# other programs
OBABEL          = "obabel"
FORCEBALANCE    = "ForceBalance"
