from Common import * 
from update import * 
from Molecule import Molecule 
from SolventBox import SolventBox 


mol = Molecule(name='perchlorate', 
               charge=-1, 
               multiplicity=2, 
               mol_type='solvent', 
               residue_name='SLV', 
               pdb='perchlorate.pdb')

mol.update()

update_mol(mol)

box = SolventBox(name='test_box')
box.add_solvent(mol)
box.update() 

