import math
import numpy as np
import pandas as pd
from rdkit import Chem

class Solute():
    r"""
    Return ideal parameters for set up solvation for given solute molecule.
    
    Parameters
    ----------
        Name : str, Required
            Solute's IUPAC name 

        SMILES String : str, Required
            Solute's SMILES string for determining its spin multiplicity via get_spin_multiplicity()
        
        Charge: int, Required
            Solute's charge
        
        File Path: str, Required
            Solute's XYZ file path for determining the ideal solvation box length using get_box_length()

    Returns
    -------
        Integer: 1 for sucessful completition; Otherwise problem occured while running.
    
    """

    def __init__(self, name, sml, charge, filePath):
        self.name = name
        self.sml = sml
        self.charge = charge
        self.file = filePath
        self.xyz = np.zeros((2,3))
        self.box_len = 0
        self.spin_multiplicity = 1
    
    def read_xyz(self):
        r"""
        Read the given XYZ file and turn it into a numpy array.

        Parameters
        ----------
            NONE

        Returns
        -------
            Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        
        """
        try:
            xyz = open(self.file, "r")
            cord = []
            for line in xyz:
                try:
                    line_list = line.split()
                    if len(line_list) == 4:
                        cord.append([float(line_list[1]), float(line_list[2]), float(line_list[3])])
                except Exception as e:
                    print(e)
                    return 0
            self.xyz = np.array(cord)
            return 1
        except:
            print(f'Error, invalid XYZ file path {self.file}.')
            return 0
    
    def distance_3D(self, x, y, z, x1, y1, z1):
        r"""
        Given x, y, z coordinates of two atoms, find the distance between them.
    
        Parameters
        ----------
            x: int, Required, x coordinate of atom1
            y: int, Required, y coordinate of atom1
            z: int, Required, z coordinate of atom1
            x1: int, Required, x coordinate of atom2
            y1: int, Required, y coordinate of atom2
            z1: int, Required, z coordinate of atom2

        Returns 
        -------
            Distance: int

        """
        return math.sqrt((x1 - x) ** 2 + (y1 - y) ** 2 + (z1 - z) ** 2)
    
    def corresponding_distance(self, column, column1, column2):
        r"""
        Find the distance between two furthest atoms in the first provided column values and their 
        corresponding coordinates to determine the distance in a 3D cartesian coordinate system.
    
        Parameters
        ----------
            col: array, Required
            col1: array, Required
            col2: array, Required

        Returns 
        -------
            Distance: int

        """
        right = np.where(column == column.max())
        left = np.where(column == column.min())

        dis = self.distance_3D(column[left], column1[left], column2[left], column[right], column1[right], \
                               column2[right])

        return dis
    
    def set_box_length(self):
        r"""
        Update the ideal solvation box length for a given solute.

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        
        """
        if not self.read_xyz(): 
            return 0
        else:
            x_cords = self.xyz[:, 0]
            y_cords = self.xyz[:, 1]
            z_cords = self.xyz[:, 2]

            distance_x = self.corresponding_distance(x_cords, y_cords, z_cords)
            distance_y = self.corresponding_distance(y_cords, x_cords, z_cords)
            distance_z = self.corresponding_distance(z_cords, y_cords, x_cords)

            dist = max([distance_x, distance_y, distance_z])
            if self.charge == 0:
                self.box_len = dist + 18
            else:
                self.box_len = dist + 26
            return 1
    
    def set_spin_multiplicity(self):
        r"""
        Update the spin mutiplicity for a given solute.

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        """
        try:
            mol = Chem.MolFromSmiles(self.sml)
            unpaired_electrons = 0
            
            for atom in mol.GetAtoms():
                try:
                    unpaired_electrons += atom.GetNumRadicalElectrons()
                except Exception as e:
                    print(e)
                    return 0
                
            total_angular_momentum = (unpaired_electrons / 2) * ((unpaired_electrons / 2) + 1)
            self.spin_multiplicity = 2 * total_angular_momentum + 1
            return 1
        except:
            print(f'Error, incorrect SMILES input {self.sml}')
            return 0

    def get_box_length(self):
        r"""
        Return the ideal solvation box length for given solute

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Solvation box length: int
        
        """
        if self.set_box_length():
            return self.box_len
        else:
            print('Unsuccessful running, please check the error messages above.')
            return
        
    def get_charge(self):
        r"""
        Return the charge of the given solute

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Charge of solute: int
        
        """
        return self.charge
    
    def get_SMILES(self):
        r"""
        Return the SMILES string of the given solute

        Parameters
        ----------
            NONE
        
        Returns
        -------
            SMILES string of solute: int
        
        """
        return self.sml

    def get_coordinates(self):
        r"""
        Return the cartesian coordinates of the given solute

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Cartesian coordinates of solute: int
        
        """
        cord = pd.DataFrame(self.xyz)
        cord.columns = ['X', 'Y', 'Z']
        return cord
    
    def get_spin_multiplicity(self):
        r"""
        Return the spin multiplicity of the given solute

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Spin multiplicity of solute: int
        
        """
        if self.set_spin_multiplicity():
            return self.spin_multiplicity
        else:
            print('Unsuccessful running, please check the error messages above.')
            return
    
    def get_methods(self):
        r"""
        Return the suggested method/methods suitable for the given solute's solvation

        Parameters
        ----------
            NONE
        
        Returns
        -------
            Suggested method/methods: list
        
        """
        if self.spin_multiplicity == 1:
            return ['resp', 'bcc']
        elif self.spin_multiplicity > 1:
            return ['resp']
        else:
            print('Error, incorrect spin multiplicity calculated.')
            return

