import os 
import sys
sys.path.append(os.getcwd().split('autosolvate/')[0] + '/autosolvate/utils/') 
from Molecule import Molecule 
from SolventBox import SolventBox 
from update import *


slu = Molecule(name='naphthalene_neutral',
               charge=0, 
               multiplicity=1, 
               mol_type='solute', 
               residue_name='SLU',
               xyz='naphthalene_neutral.xyz')

update_mol(slu)

slv = Molecule(name='dmso',
               charge=0, 
               multiplicity=1,
               mol_type='solvent', 
               residue_name='SLV',
               lib='dmso.lib', 
               frcmod='dmso.frcmod')


box = SolventBox(name='text_box')
box.add_solute(slu)
box.add_solvent(slv)
update_box(box) 
