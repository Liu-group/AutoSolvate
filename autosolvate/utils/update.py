from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np
import getopt, sys, os, subprocess, pkg_resources, glob 
from dataclasses import dataclass, field, asdict
from Common             import * 
from Molecule           import Molecule 
from SolventBox         import SolventBox 
from AntechamberDocker  import AntechamberDocker 
from TleapDocker        import TleapDocker 
from PackmolDocker      import PackmolDocker
from ParmchkDocker      import ParmchkDocker 



def update_mol(mol: object) -> None:
    r'''
    @TODO: 
        1. implement update for solute, or there is no need to update solute 
    '''
    if mol.mol_type == 'solvent': 
        update_solvent(mol)
    else: 
        print('not implemented yet')
    return 
        


def update_solvent(mol: object) -> None:
    #four cases: 

    #solvent ready in amber library
    if mol.pdb is None and mol.xyz is None: 
        if mol.box is not None and mol.frmod is not None and mol.lib is not None: 
            pass 
    

    #only has xyz ready 
    if mol.xyz is not None and mol.pdb is None:
        if mol.box is None and mol.frmod is None and mol.lib is None and mol.mol2 is None: 
            '''
            @TODO:
            1. implement xyz to pdb first 
            then starts the next iteration
            '''
            pass 

    #solvent ready in data/ folder  
    if mol.pdb is not None: 
        if mol.frcmod is not None:
            if mol.mol2 is None and mol.lib is None: 
                '''
                @TODO: 
                1. double check with Fangning if it is ok to only have frcmod file 
                ''' 
                pass

    #custome solvent only has pdb file 
    if mol.pdb is not None: 
        '''
        @TODO:
        1. custome solvent. use AntechamberDocker to convert the solvet 
        '''
        ante = AntechamberDocker() 
        ante.run(mol)
        mol.update() 
        
        parmk = ParmchkDocker() 
        parmk.run(mol) 

        


def pack_solventox(box: object) -> None:
    pack = PackmolDocker()
    pack.run(box)
    box.update()