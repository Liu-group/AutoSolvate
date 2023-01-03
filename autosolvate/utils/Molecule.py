# @TODO
# 1. modify the @EXAMPLE in the docstring 
# 2. implement charge and multiplicity determination
# 3. TleapDocker().run(mol) not fully implemented yet 
# 4. check update() method in Molecule class 

# @DISCUSSION 
# 1. should we store doc object instead of a string of the filename?
#    for mol2, frcmod, lib, prmtop, inpcrd, box
from openbabel import pybel
from openbabel import openbabel as ob
import numpy as np
import getopt, sys, os, subprocess, pkg_resources, glob, shutil
from dataclasses import dataclass, field, asdict
from Common import * 
import tools 


@dataclass 
class Molecule: 
    #constants 
    _SUPPORT_INPUT_FORMATS = ['pdb', 'xyz'] 
    
    #required arguments
    name:           str  
    charge:         int
    multiplicity:   int
    mol_type:       str

    #positional arguments 
    residue_name:   str = 'MOL'
    #files 
    pdb:            str = None
    xyz:            str = None
    mol2:           str = None
    frcmod:         str = None
    lib:            str = None 
    prmtop:         str = None
    inpcrd:         str = None
    '''
    @TODO: 
    1. this box attribute name is confusing with box as SolventBox object 
    2. we should use a better name
    '''
    box:            str = None

    folder:         str = field(default=None, init=False)
    init_folder:   bool = True

   
    def __post_init__(self) -> None:
        ''' 
        @NOTE:
        1. amber dic solvent can not use os.path.exists() to check if the frcmod exists 
        
        2. some frcmod files is named as frcmod.*  (not *.frcmod)[ask Dr.Liu]

        3. remember molecule can be initailzed many ways: 
            1. from pdb file
            2. from xyz file (not implemented yet) 
            3. from amber library. (we provide lib, frcmod and box files) 
            4. from mol2 and frcmod files (we provide these files in data/) 
 
        @TODO:
        1. custome solvent does have mol2 and frcmod files.[ask Dr.Liu]
           we need to check if mol2 and frcmod files provided in the working directory or data/ 
        '''

        if not isinstance(self.name, str): 
            raise Exception('name is not a string') 

        if self.charge not in [-1, 0, 1]: 
            raise Exception('invalid charge') 

        if self.multiplicity not in [1, 2, 3]: 
            raise Exception('invalid multiplicity')
            
        if self.mol_type not in ['solute', 'solvent']:
            raise Exception('invalid mol_type')

        if not isinstance(self.residue_name, str): 
            raise Exception('residue_name is not a string') 

        if self.pdb is not None:
            if not isinstance(self.pdb, str): 
                raise Exception('pdb is not a string') 
            if not os.path.exists(self.pdb): 
                raise Exception('pdb file does not exist') 

        if self.xyz is not None: 
            if not isinstance(self.xyz, str): 
                raise Exception('xyz is not a string') 
            if not os.path.exists(self.xyz): 
                raise Exception('xyz file does not exist') 

        if self.mol2 is not None: 
            if not isinstance(self.mol2, str): 
                raise Exception('mol2 is not a string') 
                
        if self.frcmod is not None:
            if not isinstance(self.frcmod, str): 
                raise Exception('frcmod is not a string')

        if self.lib is not None:
            if not isinstance(self.lib, str): 
                raise Exception('lib is not a string') 

        if self.prmtop is not None:
            if not isinstance(self.prmtop, str): 
                raise Exception('prmtop is not a string') 
            
        if self.inpcrd is not None: 
            if not isinstance(self.inpcrd, str): 
                raise Exception('inpcrd is not a string')

        if self.box is not None: 
            if not isinstance(self.box, str): 
                raise Exception('box is not a string') 

        if self.init_folder:
            self.set_folder()

        self.update()


    def __eq__(self, o: object) -> bool: 
        ''' 
        @TODO 
        test this method 
        ''' 
        if not isinstance(o, Molecule):
            return False 
        for key in asdict(self).keys(): 
            if getattr(self, key) != getattr(o, key): 
                return False
        return True


    def set_folder(self) -> None:
        path = WORKING_DIR + self.name + '/'
        if not os.path.exists(path): 
            os.makedirs(path) 
        self.folder = self.name + '/' 


    def update(self) -> None:
        ''' 
        @TODO
        1. can be simplified 
        ''' 
        search_range = WORKING_DIR + '/**/*'
        
        for file in glob.glob(search_range, recursive=True):
            if os.path.isfile(file):
                if file.endswith(self.name + '.pdb'):
                    self.pdb = self.name + '.pdb'
                    #copy file to mol/ folder 
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                
                elif file.endswith(self.name + '.xyz'):
                    self.xyz = self.name + '.xyz'
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/') 
                
                elif file.endswith(self.name + '.mol2'):
                    self.mol2 = self.name + '.mol2'
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')
                
                elif file.endswith(self.name + '.frcmod'):
                    self.frcmod = self.name + '.frcmod' 
                    if self.name+'/' not in file: 
                        shutil.copy(file, WORKING_DIR + self.name + '/') 

                elif file.endswith(self.name + '.lib'):
                    self.lib = self.name + '.lib'
                    if self.name+'/' not in file:
                        shutil.copy(file, WORKING_DIR + self.name + '/')
                


#AMBER SOLVENTS 
AMBER_WATER               = Molecule(name='water', charge=0, multiplicity=1, 
                                     mol_type='solvent', 
                                     box='TIP3PBOX', 
                                     init_folder=False 
                            )

AMBER_METHANOL            = Molecule(name='methanol', charge=0, multiplicity=1, 
                                     mol_type='solvent', 
                                     lib='solvents.lib', frcmod='frcmod.meoh', box='MEOHBOX',
                                     init_folder=False
                            )

AMBER_CHLOROFORM          = Molecule(name='chloroform', charge=0, multiplicity=1, 
                                     mol_type='solvent',  
                                     lib='solvents.lib', frcmod='frcmod.chcl3', box='CHCL3BOX', 
                                     init_folder=False
                            )
                            
AMBER_NMA                 = Molecule(name='nma', charge=0, multiplicity=1, 
                                     mol_type='solvent',
                                     lib='solvents.lib', frcmod='frcmod.nma', box='NMABOX', 
                                     init_folder=False 
                            )

AMBER_SOLVENT_LIST       =  [AMBER_WATER, AMBER_METHANOL, AMBER_CHLOROFORM, AMBER_NMA] 