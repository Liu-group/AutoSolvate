"""
resp_factory.py
Automated workflow to to run resp fitting with different quantum chemistry packages

Handles the primary functions
"""
from autosolvate.resp_classes.resp_gaussian import RespGaussian
from autosolvate.resp_classes.resptools.resp_gamess import RespGAMESS
from autosolvate.resp_classes.resp_orca import RespORCA
from  autosolvate.resp_classes.resp_orca import get_crds_from_xyz

def resp_factory(**kwargs):
    """ 
    factory class for calculating diagnostics

    Parameters
    ----------
    required keyword arguments:
         pdbfile : str
             file name of the pdb file for RESP fitting
         
         qm_program: str
             Quantum chemistry pacakge used to calculat RESP
             Currently support: ['gamess','gaussian','orca']

    Returns
    -------
    cls_instance: instance of the class that performs RESP fitting with the request quantum chemistry program
    """
    cls_dict = dict(gaussian=RespGaussian, gamess=RespGAMESS, orca=RespORCA)

    qm_program = kwargs["qm_program"] if "qm_program" in kwargs.keys() else "gamess"

    if qm_program not in cls_dict.keys():
        raise Exception("The request quantum chemistry program, {}, is not supported.".format(qm_program))

    if "pdbfile" not in kwargs.keys():
        raise KeyError("Error: funciton resp_factory requeests a pdbfile file in the keyword argument for resp charge fitting.")
    

    cls = cls_dict[qm_program]
    cls_instance = cls(**kwargs)
    return cls_instance
