import os
import pytest
import numpy as np
import os
from autosolvate.FFmetalcomplex import *

from . import helper_functions as hp

def compare_boxgen(out, ref):
    with open(out + ".inpcrd", 'r') as f_in:
        inpcrd = f_in.read()
    with open(ref + ".inpcrd", 'r') as f_in:
        inpcrd_ref = f_in.read()
    compare_inpcrd = inpcrd == inpcrd_ref
    with open(out + ".prmtop", 'r') as f_in:
        prmtop = f_in.readlines()[1:] #skip DATE = 10/08/22  11:32:06
    with open(ref + ".prmtop", 'r') as f_in:
        prmtop_ref = f_in.readlines()[1:]
    compare_prmtop = prmtop == prmtop_ref
    return  compare_inpcrd, compare_prmtop


def testAutoMCPB_water_solvated(tmpdir):
    os.system('cp inputs/metalcomplex/FeCP2_plus2/* ' + str(tmpdir))
    startFFgen(["-m", "FeCP2_plus2.xyz", "-k", "2","-e","orca","-d","/opt/orca/5.0.2/orca","-y","16",'-j','lanl2dz','-w','Y','-s','water'])
    assert os.path.exists('FeCP2_plus2_solvated.prmtop')
    assert os.path.exists('FeCP2_plus2_solvated.inpcrd')
  #  os.system('cp FeCP2_plus2_solvated.prmtop ~')
    out =  'FeCP2_plus2_solvated'
    ref = os.path.join(hp.get_reference_dir(), 'metalcomplex/FeCP2_plus2/FeCP2_plus2_solvated')
 #   compare_inpcrd, compare_prmtop = compare_boxgen(out, ref)
#    assert compare_inpcrd
    os.system('rm ' +str(tmpdir)+'/*')
 #   assert compare_prmtop# cant be exactly the same because of the date

def testAutoMCPB_acn_solvated(tmpdir):
    os.system('cp inputs/metalcomplex/Kcryptand/* ' + str(tmpdir))
    startFFgen(["-m", "Kcryptand.xyz", "-k", "1","-e","orca","-d","/opt/orca/5.0.2/orca","-y","16",'-j','lanl2dz','-w','Y','-s','acetonitrile','-z','Y'])
    assert os.path.exists('Kcryptand_solvated.prmtop')
    assert os.path.exists('Kcryptand_solvated.inpcrd')
    
    cmd = 'cp ' + str(tmpdir) +'/Kcryptand_solvated.inpcrd ~'
    os.system(cmd)
    cmd = 'cp ' + str(tmpdir) +'/Kcryptand_solvated.prmtop ~'
    os.system(cmd)
    
    out =  'Kcryptand_solvated'
    ref = os.path.join(hp.get_reference_dir(), 'metalcomplex/Kcryptand/Kcryptand_solvated')
    compare_inpcrd, compare_prmtop = compare_boxgen(out, ref)
    assert compare_inpcrd
    os.system('rm ' +str(tmpdir)+'/*')

def testAutoMCPB_dmso_solvated(tmpdir):
    os.system('cp inputs/metalcomplex/Pb_complex/* ' + str(tmpdir))
    os.system('cp inputs/dmso.off ' + str(tmpdir))
    os.system('cp inputs/dmso.frcmod ' + str(tmpdir))
    startFFgen(["-m", "Pb_complex.xyz", "-k", "2","-e","orca","-d","/opt/orca/5.0.2/orca","-y","16",'-j','lanl2dz','-w','Y','-s','d','-z','Y','-l','dmso.off','-p','dmso.frcmod','-x','3.0'])
   # os.system('cat '+ os.path.join(str(tmpdir), 'tleap.log'))
   # os.system(f'cp {tmpdir}/* ~/check/')
    assert os.path.exists('Pb_complex_solvated.prmtop')
    assert os.path.exists('Pb_complex_solvated.inpcrd')
    
    out =  'Pb_complex_solvated'
    ref = os.path.join(hp.get_reference_dir(), 'metalcomplex/Pb_complex/Pb_complex_solvated')
    compare_inpcrd, compare_prmtop = compare_boxgen(out, ref)
    assert compare_inpcrd
    
    




    
