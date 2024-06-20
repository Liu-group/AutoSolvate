import pubchempy as pcp
from openbabel import pybel
import os

class PubChemAPI():
    r"""
    Solute molecule construction using PubChem API.
    
    Parameters
    ----------
        Name : str, Required
            Given solute's IUPAC name, solute's information is obtained from PubChem API. 
        
            Simplified molecular-input line-entry system (SMILES) form and charge can
            be retreived using get_info().

        FilePath : str, Optional, default: local directory
            Solute's XYZ file generated using openbabel is saved to given file path.

    Returns
    -------
        Integer: 1 for sucessful completition; Otherwise problem occured while running.
    
    """

    def __init__(self, name, filePath = os.getcwd()):
        self.name = name
        self.cID = ''
        self.mol = None
        self.sml = 'NONE'
        self.charge = 0
        self.ff = 'mmff94'
        self.s = 50
        self.filePath = filePath
        self.export = os.path.join(self.filePath, 'tempFile')

    def set_up_mol(self):
        r"""
        Set up molecule information using PubChem API.
        
        Parameters
        ----------
            NONE

        Returns
        -------
            Integer: 1 for sucessful completition; Otherwise problem occured while running.
        
        """
        try:
            self.cID = pcp.get_cids(self.name)
            try:
                self.mol = pcp.Compound.from_cid(self.cID)
                return 1
            except:
                print(f'Error, PubChem API could not access information about given compound ID {self.cID} \
                      from the molecule database')
                return 0
        except:
            print(f'Error, {self.name} corresponding compound ID cannot be obtained from PubChem API.')
            return 0


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
            self.sml = self.mol.canonical_smiles
            return 1
        except:
            print(f'Error, solute molecule is not correctly set up, make sure set_up_mol() is called before.')
            return 0
    
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
            self.charge =  self.mol.charge
            return 1  
        except:
            print(f'Error, solute molecule is not correctly set up, make sure set_up_mol() is called \
                  before.')
            return 0

    def set_force_field(self, forcefield = 'mmff94'):
        r"""
        Set force field used in 3-D coordinates construction.
        
        Parameters
        ----------
            Force field: str, Optional, default: 'mmff94'

            Available force field: 'mmff94', 'uff', 'ghemical'
        
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
            Steps: int, steps used in 3-D coordniates construction.
        
        """
        return self.s
    
    def set_path(self, path = os.getcwd()):
        r"""
        Set file path to save solute's xyz file.
        
        Parameters
        ----------
            path: str, Optional, default: local directory
        
        Returns
        -------
            NONE
        
        """
        self.filePath = path
    
    def get_path(self):
        r"""
        Return file path where solute's xyz file is saved.
        
        Parameters
        ----------
            NONE
        
        Returns
        -------
            File path: str
        
        """
        return self.filePath
    
    def get_XYZ(self):
        r"""
        Save XYZ file of the solute molecule.

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        
        """
        try:
            compound_infos = pcp.get_compounds(self.name, 'name', record_type='3d')
            try:
                coord_info = compound_infos[0]
                atoms = coord_info.atoms
                num_atoms = str(len(atoms)) + '\n'
                coordinates_3d = [(atom.element, atom.x, atom.y, atom.z) for atom in atoms]

                with open(self.export, 'w') as file:
                    file.write(num_atoms)
                    for line in coordinates_3d:
                        coord = '\n '
                        for elem in line:
                            if type(elem) == str:
                                coord += ' '
                                coord += elem
                            else:
                                coord += ' '
                                coord += str(elem)
                        file.write(coord)
                print(f'XYZ file successfully saved as {self.export}')
                return 1
            except:
                print('3D coordinate is unavailable from PubChem, try upload a xyz file')
                return 0
        except:
            print(f'Error, invalid IUPAC name: {self.name}')
            return 0
    
    def get_info(self):
        r"""
        Provide information for given solute

        Parameters
        ----------
            NONE
        
        Returns
        -------
            (IUPAC Name: str, SMILES: str, Charge: int, XYZ File Path: str) : Tuple
        
        """
        if self.set_up_mol() and self.set_SMILES() and self.set_charge() and self.get_XYZ():
            return (self.name, self.sml, self.charge, self.export)
        else:
            print('Unsuccessful running, please check the error messages above.')
            return
