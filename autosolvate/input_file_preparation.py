import pubchempy as pcp
from openbabel import pybel

def get_molecule(name):
    r"""
    Turn IUPAC name of a molecule into SMILES format
   
    Parameters
    ----------
    name: str, Required, default: None
          IUPAC name of a molecule

    Returns 
    -------
    SMILES string representation of provided molecule

    """
    try:
        sml = pcp.get_properties('CanonicalSMILES', name,'name')[0]['CanonicalSMILES']
        return sml
    except:
        print(f'Cannot find {name} in PubChem data base, please provide IUPAC name')
        return
    
def smile_to_xyz(sml, ff = 'mmff94', s = 50):
    r"""
    Turn SMILES string to XYZ file
   
    Parameters
    ----------
    sml: str, Required, default: None
         SMILES string
    ff:  str, Optional, default: mmff94
         Use mmff94 or uff or ghemical
    s:   int, Optional, default: 50
         Steps greater or equal to 50

    Returns
    -------
    2D structure of the SMILES molecule

    """
    try:
        mol = pybel.readstring("smi", sml)
        try:
            mol.make3D(forcefield=ff, steps=s)
            mol.write('xyz', "mol.xyz", overwrite=True)
            return mol
        except:
            print(f'Could not convert {sml} into 3D structure with forcefield {ff} with steps {s}.\nTry a different parameter setup.')
            return
    except:
        print('Invalid SMILES input, please try again.')
        return