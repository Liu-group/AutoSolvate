import os
import shutil

import pytest

import autosolvate
from autosolvate.ion_forcefield import IonForceFieldBuilder


_REQUIRED = ("obabel", "tleap")
pytestmark = pytest.mark.skipif(
    any(shutil.which(exe) is None for exe in _REQUIRED),
    reason="AmberTools and Open Babel are required for ion integration tests",
)


def _assert_amber_outputs(prefix, element, charge):
    parmed = pytest.importorskip("parmed")
    for suffix in (".pdb", ".prmtop", ".inpcrd"):
        assert os.path.isfile(prefix + suffix)
    topology = parmed.load_file(prefix + ".prmtop", prefix + ".inpcrd")
    assert any(atom.element_name.upper() == element.upper() or atom.name.upper() == element.upper() for atom in topology.atoms)
    counterion = "Cl" if charge > 0 else "Na"
    assert sum(atom.element_name == counterion for atom in topology.atoms) == abs(charge)


def test_potassium_water():
    autosolvate.startboxgen([
        "-m", "inputs/K.xyz", "-s", "water", "-c", "1", "-o", "K_solvated",
    ])
    _assert_amber_outputs("K_solvated", "K", 1)
    frcmod = open("inputs/K.frcmod").read()
    assert "MASS" in frcmod and "NONB" in frcmod


def test_sulfur_acetonitrile():
    if shutil.which("packmol") is None:
        pytest.skip("Packmol is required for acetonitrile")
    autosolvate.startboxgen([
        "-m", "inputs/S.xyz", "-s", "acetonitrile", "-c", "-2", "-o", "S_solvated",
    ])
    _assert_amber_outputs("S_solvated", "S", -2)


def test_iron_custom_dmso():
    autosolvate.startboxgen([
        "-m", "inputs/Fe.xyz", "-s", "d", "-c", "2", "-o", "Fe_dmso_solvated",
        "-l", "inputs/dmso.off", "-p", "inputs/dmso.frcmod",
    ])
    _assert_amber_outputs("Fe_dmso_solvated", "Fe", 2)
