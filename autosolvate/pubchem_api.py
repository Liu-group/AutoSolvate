import pubchempy as pcp

class solute_from_PubChem():
    r"""
    Solute molecule construction using PubChem API.
    
    Parameters
    ----------
    Name : str, Required
        Given solute's IUPAC name, solute's information is obtained from PubChem API. 
        1) Simplified molecular-input line-entry system (SMILES) form can be accessed 
        by get_SMILES(). 2) Solute's charge is retreived using get_charge().
    
    Returns
    -------
    Integer: 1 for sucessful completition; Otherwise problem occured while running.
    
    """


    def __init__(self, name = None):
        self.name = name
        self.sml = 'NONE'
        self.ff = 'mmff94'
        self.s = 50
        self.charge = 0

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
    
    def get_force_field(self):
        r"""
        Return the force field used in 3D coordinate generation.
        
        Parameters
        ----------
            NONE

        Returns
        -------
            Force field: str, force field used in the 3D coordinates 
            construction
        
        """
        return self.ff
    
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
    
    def get_steps(self):
        r"""
        Return steps used in 3-D coordinates construction. 
        
        Parameters
        ----------
            NONE
        
        Returns
        -------
            Steps: int, steps used in 3-D coordniates 
            construction.
        
        """
        return self.s

    def set_SMILES(self):
        r"""
        Update SMILES form of the solute molecule using PubChem API.

        Parameters
        ----------
            NONE
        
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
            NONE
        
        Returns
        -------
            SMILES: str
        
        """
        return self.sml

    def set_charge(self):
        r"""
        Update charge of the molecule using PubChem API
    
        Parameters
        ----------
            NONE

        Returns 
        -------
            Integer: 1 for sucessful completition; 0 for imcompletition.

        """
        try:
            cID = pcp.get_cids(self.name)
            molecule = pcp.Compound.from_cid(cID)
            self.charge =  molecule.charge
            return 1  
        except:
            print(f'Cannot find {self.name} in PubChem data base, please provide IUPAC name')
            return 0

    def get_charge(self):
        r"""
        Return solute's charge
    
        Parameters
        ----------
        None

        Returns 
        -------
        charge: int

        """
        return self.charge
