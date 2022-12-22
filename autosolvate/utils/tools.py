r''' 
    @TODO 
    1. check '-> callable:' is used correctly 
'''
from Common import * 
import os   
from openbabel import openbabel as ob


def extract_basename_name_extension(inputfile: str) -> tuple:
    r'''
    Get basename, name and extension of a file

    @Example: 
    >>> extract_basename_name_extension('test.pdb') 
    ('test.pdb', 'test', 'pdb')
    '''
    basename        = os.path.basename(inputfile) 
    name , ext      = os.path.splitext(basename)
    ext             = ext[1:]
    return basename, name, ext


def convert_xyz_to_pdb(inputfile: str) -> str: 
    r'''
    Convert xyz file to pdb file using openbabel
    '''
    basename, name, ext = extract_basename_name_extension(inputfile) 
    if ext == 'xyz':
        outputfile = name + '.pdb'
        os.system('ob -ixyz {} -opdb -O {}'.format(inputfile, outputfile))
        return outputfile
    else:
        raise Exception('Input file format not supported in function convert_xyz_to_pdb()')


def srun() -> callable:
    r'''
    @Example:
    >>> @srun()
    ... def test():
    ...     return 'cmd'
    >>> test()
    'srun -n 1 cmd'
    '''
    def wrap(func): 
        def wrapped_func(*args, **kwargs): 
            if USE_SRUN:  
                value = 'srun -n 1 ' + func(*args, **kwargs) 
            else:
                value = func(*args, **kwargs)  
            return value
        return wrapped_func 
    return wrap

            

if __name__ == '__main__': 
    import doctest

    global USE_SRUN, DRY_RUN
    USE_SRUN = True
    
    doctest.testmod() 