#@TODO: 
#       is a better way to implement this? : central control of the whole project

from Molecule import Molecule
import os 


#global variable for the whole project 
DRY_RUN         = False
USE_SRUN        = False


WORKING_DIR     = os.getcwd() + '/'



AMBER_WATER               = Molecule(name='water', charge=0, multiplicity=1, 
                                     mol_type='solvent', 
                                     box='TIP3PBOX'
                            )

AMBER_METHANOL            = Molecule(name='methanol', charge=0, multiplicity=1, 
                                     mol_type='solvent', 
                                     lib='solvents.lib', frcmod='frcmod.meoh', box='MEOHBOX'
                            )

AMBER_CHLOROFORM          = Molecule(name='chloroform', charge=0, multiplicity=1, 
                                     mol_type='solvent',  
                                     lib='solvents.lib', frcmod='frcmod.chcl3', box='CHCL3BOX'
                            )
                            
AMBER_NMA                 = Molecule(name='nma', charge=0, multiplicity=1, 
                                     mol_type='solvent',
                                     lib='solvents.lib', frcmod='frcmod.nma', box='NMABOX'
                            )

AMBER_SOLVENT_LIST       =  [AMBER_WATER, AMBER_METHANOL, AMBER_CHLOROFORM, AMBER_NMA] 