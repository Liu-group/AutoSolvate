## @file __main__.py
#  Gateway script to rest of program
#
# !/usr/bin/env python

from autosolvate.GUI.tk_autosolvate import *
from autosolvate.autosolvate import *
from autosolvate.generatetrajs import *
from autosolvate.clustergen import *
from autosolvate.FFmetalcomplex import *
from autosolvate.multicomponent import *

## Main function
#  @param args Argument namespace
def main(args=None):
    if args is None:
        args = sys.argv[1:]
    ### run with gui ###
    if len(args) == 0:
        print('AutoSolvate is starting with graphical user interface!')
        ### create main application
        window = Tk()
        my_gui = autosolvateGUI(window)
        window.mainloop()
        cleanUp()
    elif args[0] == '-h' or args[0] == '--help':
        print('Usage: autosolvate [OPTIONS]')
        print('  [NO OPTION]                launches GUI')
        print('  boxgen [OPTIONS]           generate initial structure')
        print('  boxgen_metal[OPTIONS]      generate initial structure for organometallic compounds')
        print('  mdrun [OPTIONS]            automated QM/MM trajectory generatio')
        print('  clustergen [OPTIONS]       extract microsolvated clusters ')
        print('  -h, --help                 short usage description')
        print()
        print('All options described at autosolvate.readthedocs.io')
        exit()
 
    elif args[0] == 'boxgen':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to generate solvent box and force field parameters.')
        startboxgen(args[1:])
    elif args[0] == 'boxgen_metal':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to generate solvent box and force field parameters for organometallic compounds.')
        startFFgen(args[1:])
    elif args[0] == 'mdrun':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to automatically run MD simulations of solvated structure.')
        startmd(args[1:])
    elif args[0] == 'clustergen':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to extract microsolvated clusters from MD trajectories with solvent box.')
        startclustergen(args[1:])
    elif args[0] == 'boxgen_multicomponent':
        print('AutoSolvate is starting in command line mode!')
        print('Running the module to generate solvent box and force field parameters for multicomponent systems.')
        startmulticomponent(args[1:]) 
    else:
        print('Invalid syntax for AutoSolvate command line interface.')
        print('Please run \'autosolvate -h\' to check out the basic usage.')
        exit()

if __name__ == '__main__':
    main()
