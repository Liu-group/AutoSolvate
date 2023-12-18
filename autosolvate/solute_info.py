import math
from openbabel import pybel

def distance_3d(x, y, z, x1, y1, z1):
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


class solute_info():
    r"""
    Access solute molecule information.
    
    Parameters
    ----------
    Sml: str, Required
        Given solute's SMILES string. Its XYZ file is generated using set_XYZ() 
        function. The XYZ file can be used to obtain ideal solvation box length
        via get_box_length()
    
    Returns
    -------
    Integer: 1 for sucessful completition; Otherwise problem occured while running.
    
    """

    def __init__(self, sml = None):
        self.sml = sml
        self.filePath = ''
        self.export = self.filePath + "mol.xyz"
        self.box_len = 0
    
    def set_path(self, path):
        r"""
        Set file path to save solute's xyz file.
        
        Parameters
        ----------
        path: str, Optional, default: ''
        
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
        None
        
        Returns
        -------
        Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        
        """
        try:
            mol = pybel.readstring("smi", self.sml)
            try:
                mol.make3D(forcefield=self.ff, steps=self.s)
                mol.write('xyz', self.export, overwrite=True)
                print(f'XYZ file successfully saved to {self.filePath}')
                return 1
            except:
                print(f'Error, could not convert {self.sml} into 3D structure with forcefield {self.ff} with steps {self.s}.\nTry a different parameter setup.')
                return 0
        except:
            print('Error, invalid SMILES input, please try again.')
            return 2
    
    def set_box_length(self, charge):
        r"""
        Update the ideal solvation box length for a given solute.

        Parameters
        ----------
        charge: int, Required
            Charged molecule will have more solvation box length compare to neutral solute.
        
        Returns
        -------
        Integer: 1 for sucessful completition; Otherwise, it is for imcompletition.
        
        """
        try:
            file = self.export
            xyz = open(file, "r")
        
            x_list = []
            y_list = []
            z_list = []

            for line in xyz:
                line_list = line.split()
                if len(line_list) == 4:
                    x_list.append(float(line_list[1]))
                    y_list.append(float(line_list[2]))
                    z_list.append(float(line_list[3]))
            
            try:
                xmin = min(x_list)
                ymin = min(y_list)
                zmin = min(z_list)

                xmax = max(x_list)
                ymax = max(y_list)
                zmax = max(z_list)
                
                left_ix = x_list.index(xmin)
                right_ix = x_list.index(xmax)

                left_iy = y_list.index(ymin)
                right_iy = y_list.index(ymax)

                left_iz = z_list.index(zmin)
                right_iz = z_list.index(zmax)

                distance_x = distance_3d(x_list[left_ix], y_list[left_ix], z_list[left_ix], x_list[right_ix], y_list[right_ix], z_list[right_ix])
                distance_y = distance_3d(x_list[left_iy], y_list[left_iy], z_list[left_iy], x_list[right_iy], y_list[right_iy], z_list[right_iy])
                distance_z = distance_3d(x_list[left_iz], y_list[left_iz], z_list[left_iz], x_list[right_iz], y_list[right_iz], z_list[right_iz])

                dist = max([distance_x, distance_y, distance_z])

                if charge == 0:
                    self.box_len = dist + 18
                else:
                    self.box_len = dist + 26
                return 1
            except:
                print('Error, incorrect XYZ file format, supported format has form:')
                print('ATOM          X        Y       Z')
                return 0
        except:
            print(f'Error, file cannot be read, check if file exist at {self.export}')
            return 0

    def get_box_length(self):
        r"""
        Return the ideal solvation box length for given solute

        Parameters
        ----------
        None
        
        Returns
        -------
        Solvation box length: int
        
        """
        return self.box_len