import pubchempy as pcp
from openbabel import pybel

class solute_prep():
    r"""
    Solute molecule construction.
    
    Parameters
    ----------
    Name : str, Required
        Given IUPAC name, simplified molecular-input line-entry system (SMILES) form 
        can be accessed by first set_SMILES(), and then get_SMILES(). Molecule's XYZ 
        file is generated using get_XYZ() function. 
    
    Returns
    -------
    Integer: 1 for sucessful completition; Otherwise problem occured while running.
    
    """


    def __init__(self, name = None):
        self.name = name
        self.sml = 'NONE'
        self.ff = 'mmff94'
        self.s = 50
        self.path = ''

    def set_force_field(self, forcefield = 'mmff94'):
        r"""
        Set force field used in 3-D coordinates construction.
        
        Parameters
        ----------
        forcefield: str, Optional, default: 'mmff94'

        Available forcefield: 'mmff94', 'uff', 'ghemical'
        
        Returns
        -------
            NONE
        
        """
        self.ff = forcefield
    
    def set_steps(self, steps = '50'):
        r"""
        Set steps used in 3-D coordinates construction.
        
        Parameters
        ----------
        steps: int, Optional, default: '50'
        
        Returns
        -------
            NONE
        
        """
        self.s = steps
    
    def set_path(self, path = ''):
        r"""
        Set file path to save solute's xyz file.
        
        Parameters
        ----------
        path: str, Optional, default: ''
        
        Returns
        -------
            NONE
        
        """
        self.path = path

    def set_SMILES(self):
        r"""
        Update SMILES form of the solute molecule.

        Parameters
        ----------
        None
        
        Returns
        -------
        Integer: 1 for sucessful completition; 0 for imcompletition.
        
        """
        try:
            self.sml = pcp.get_properties('CanonicalSMILES', self.name,'name')[0]['CanonicalSMILES']
            return 1
        except:
            if self.name == 'NONE':
                print('Error, solute IUPAC name is not provided, please provide IUPAC name')
            else:
                print(f'Error, cannot find {self.name} in PubChem data base, please provide correct IUPAC name')
            return 0
        
    def get_SMILES(self):
        r"""
        Retreive SMILES form of the solute molecule.

        Parameters
        ----------
        None
        
        Returns
        -------
        SMILES: str
        
        """
        return self.sml

    def get_XYZ(self):
        r"""
        Save XYZ file of the solute molecule.

        Parameters
        ----------
        None
        
        Returns
        -------
        Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        
        """
        try:
            mol = pybel.readstring("smi", self.sml)
            try:
                mol.make3D(forcefield=self.ff, steps=self.s)
                mol.write('xyz', self.path + "mol.xyz", overwrite=True)
                print(f'XYZ file successfully saved to {self.path}')
                return 1
            except:
                print(f'Error, could not convert {self.sml} into 3D structure with forcefield {self.ff} with steps {self.s}.\nTry a different parameter setup.')
                return 0
        except:
            print('Error, invalid SMILES input, please try again.')
            return 2