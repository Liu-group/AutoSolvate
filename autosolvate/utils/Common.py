#@TODO: 
#       is a better way to implement this? : central control of the whole project

from Molecule import Molecule 


#global variable for the whole project 
DRY_RUN         = False
USE_SRUN        = False



WORKING_DIR     = os.getcwd() + '/'



#@TODO, move custom solvents to a Tleap class later 

AMBER_SOLVENT_DIC   = { 'water':     [' ',  'TIP3PBOX '],
                        'methanol':  ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
                        'chloroform':['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
                        'nma':       ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}


water               =  Molecule(inputfile='none', name='water', charge=0, multiplicity=1, mol_type='amber_solvent')
water.frcmod        =  'loadOff solvents.lib\n loadamberparams frcmod.meoh\n'

