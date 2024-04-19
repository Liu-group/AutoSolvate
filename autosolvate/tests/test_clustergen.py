import autosolvate
import pytest
from pkg_resources import resource_filename, Requirement
from . import helper_functions as hp

import numpy as np

def read_xyz(file):
    # return elements, coordinates
    with open(file, 'r') as f_in:
        lines = f_in.readlines()
    n_atoms = int(lines[0])
    elements = []
    coordinates = np.zeros((n_atoms, 3))
    for i in range(2, n_atoms+2):
        lines[i] = lines[i].strip()
        element, x, y, z = lines[i].split()
        elements.append(element)
        coordinates[i-2] = np.array([float(x), float(y), float(z)])
    return elements, coordinates

def compare_cluster(out, ref):
    elements_out, coordinates_out = read_xyz(out)
    elements_ref, coordinates_ref = read_xyz(ref)
    coordinates_out -= np.average(coordinates_out, axis=0)
    coordinates_ref -= np.average(coordinates_ref, axis=0)
    return np.allclose(coordinates_out, coordinates_ref, atol=1e-3) and (elements_out == elements_ref)


def test_clustergen_aspherical(tmpdir):
    autosolvate.startclustergen(["-f", hp.get_input_dir("water_solvated.prmtop"), "-t", hp.get_input_dir("water_solvated-qmmmnvt.netcdf")])
    assert compare_cluster(hp.get_input_dir("water_solvated-cutoutn-0.xyz"), hp.get_reference_dir("water_solvated-cutoutn-0.xyz")) 

def test_clustergen_spherical(tmpdir):
    autosolvate.startclustergen(["-f", hp.get_input_dir("water_solvated.prmtop"), "-t", hp.get_input_dir("water_solvated-qmmmnvt.netcdf"), "-a", "0", "-i", "10", "-s", "6", "-p"])
    assert compare_cluster(hp.get_input_dir("water_solvated-cutoutn-90.xyz"), hp.get_reference_dir("water_solvated-cutoutn-90.xyz"))
    
