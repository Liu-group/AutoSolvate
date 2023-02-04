import os 
import sys
sys.path.append(os.getcwd().split('autosolvate/')[0] + '/autosolvate/utils/') 
from Molecule import Molecule 
from Molecule import * 
from SolventBox import SolventBox 
from update import *


mol1 = Molecule(name='naphthalene_neutral', 
                charge=0, 
                multiplicity=1, 
                mol_type='solute', 
                residue_name='SLU', 
                xyz='naphthalene_neutral.xyz'
)

mol2 = Molecule(name='naphthalene_radical', 
                charge=0, 
                multiplicity=1, 
                mol_type='solute', 
                residue_name='SLU', 
                xyz='naphthalene_radical.xyz'
)


update_mol(mol1)
water = AMBER_WATER 
box = SolventBox(name='box_neutral')
box.add_solute(mol1)
box.add_solvent(water)
update_box(box)



update_mol(mol2)
water = AMBER_WATER 
box = SolventBox(name='box_radical')
box.add_solute(mol2)
box.add_solvent(water)
update_box(box)
                