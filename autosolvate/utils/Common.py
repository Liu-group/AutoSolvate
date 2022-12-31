#@TODO: 
#       is a better way to implement this? : central control of the whole project
#       1. DRY_RUN: is not working 
import os 


#global variable for the whole project 
DRY_RUN         = False
USE_SRUN        = False

WORKING_DIR     = os.getcwd() + '/'