from Common import * 
from update import * 
from Molecule import Molecule 




mol = Molecule(name='perchlorate', 
               charge=-1, 
               multiplicity=2, 
               mol_type='solvent', 
               residue_name='SLV', 
               pdb='perchlorate.pdb')

mol.update()

update_mol(mol) 