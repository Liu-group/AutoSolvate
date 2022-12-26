#@TODO: 
#       is a better way to implement this? : central control of the whole project

from Molecule import Molecule
import os 


#global variable for the whole project 
DRY_RUN         = False
USE_SRUN        = False


WORKING_DIR     = os.getcwd() + '/'



# AMBER_SOLVENT_DIC   = { 'water':     [' ',  'TIP3PBOX '],
#                         'methanol':  ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
#                         'chloroform':['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
#                         'nma':       ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}


AMBER_WATER               =  Molecule(inputfile='none', name='water', charge=0, multiplicity=1, mol_type='amber_solvent')
AMBER_WATER.lib           =  ''
AMBER_WATER.frcmod        =  ''
AMBER_WATER.box           =  'TIP3PBOX'

AMBER_METHANOL            =  Molecule(inputfile='none', name='methanol', charge=0, multiplicity=1, mol_type='amber_solvent') 
AMBER_METHANOL.lib        =  'solvents.lib' 
AMBER_METHANOL.frcmod     =  'frcmod.meoh'
AMBER_METHANOL.box        =  'MEOHBOX' 

AMBER_CHLOROFORM          =  Molecule(inputfile='none', name='chloroform', charge=0, multiplicity=1, mol_type='amber_solvent') 
AMBER_CHLOROFORM.lib      =  'solvents.lib' 
AMBER_CHLOROFORM.frcmod   =  'frcmod.chcl3' 
AMBER_CHLOROFORM.box      =  'CHCL3BOX' 

AMBER_NMA                 =  Molecule(inputfile='none', name='nma', charge=0, multiplicity=1, mol_type='amber_solvent') 
AMBER_NMA.lib             =  'solvents.lib' 
AMBER_NMA.frcmod          =  'frcmod.nma' 
AMBER_NMA.box             =  'NMABOX' 


AMBER_SOLVENT_LIST       =  [AMBER_WATER, AMBER_METHANOL, AMBER_CHLOROFORM, AMBER_NMA] 