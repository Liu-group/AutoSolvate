import autosolvate
import pytest
from pkg_resources import resource_filename, Requirement
from . import helper_functions as hp

def compare_cluster(out, ref):
    with open(out, 'r') as f_in:
        cluster_test = f_in.read()[2:]
    with open(ref, 'r') as f_in:
        cluster_ref = f_in.read()[2:]
    return cluster_test == cluster_ref

def test_clustergen_aspherical(tmpdir):
    autosolvate.startclustergen(["-f", hp.get_input_dir("water_solvated.prmtop"), "-t", hp.get_input_dir("water_solvated-qmmmnvt.netcdf")])
    assert compare_cluster(hp.get_input_dir("water_solvated-cutoutn-0.xyz"), hp.get_reference_dir("water_solvated-cutoutn-0.xyz")) 

def test_clustergen_spherical(tmpdir):
    autosolvate.startclustergen(["-f", hp.get_input_dir("water_solvated.prmtop"), "-t", hp.get_input_dir("water_solvated-qmmmnvt.netcdf"), "-a", "0", "-i", "10", "-s", "6", "-p"])
    assert compare_cluster(hp.get_input_dir("water_solvated-cutoutn-90.xyz"), hp.get_reference_dir("water_solvated-cutoutn-90.xyz"))
    
