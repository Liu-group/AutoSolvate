import os 
import sys
sys.path.append(os.getcwd().split('autosolvate/')[0] + '/autosolvate/utils/') 
from Molecule import Molecule 
from Molecule import * 
from SolventBox import SolventBox 
from update import *

slu = Molecule(name='naphthalene_neutral',
                charge=0,
                multiplicity=1,
                mol_type='solute',
                residue_name='SLU',
                xyz='naphthalene_neutral.xyz')
update_mol(slu) 

slv = Molecule(name='acetonitrile', 
                charge=0, 
                multiplicity=1, 
                mol_type='solvent', 
                residue_name='SLV', 
                pdb = 'acetonitrile.pdb', 
                prep = 'acetonitrile.prep')
update_mol(slv)

box = SolventBox(name='test_box') 
box.add_solute(slu) 
box.add_solvent(slv) 

box.duplicate_solute_num  = 1 
box.duplicate_solvent_num = 1680
box.closeness = 2.0
box.cubesize  = 56

pack_box(box)
tleap_box(box)
