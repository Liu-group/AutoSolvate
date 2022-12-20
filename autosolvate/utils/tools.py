from Common import * 
import os   


def get_basename_name_ext(inputfile):
    r'''
    Get basename, name and extension of a file

    @Example: 
    >>> get_basename_name_ext('test.pdb') 
    ('test.pdb', 'test', 'pdb')
    '''
    basename        = os.path.basename(inputfile) 
    name , ext      = os.path.splitext(basename)
    ext             = ext[1:]
    return basename, name, ext
    

def srun():
    r'''
    @Example:
    >>> @srun()
    ... def test():
    ...     return 'FirePunch'
    >>> test()
    'srun -n 1 FirePunch'
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