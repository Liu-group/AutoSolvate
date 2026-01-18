import os
from pathlib import Path
import pytest
import numpy as np 

from .helper_functions import get_input_dir

from autosolvate.resp_classes.resp_orca import RespORCA
from autosolvate.multicomponent import startmulticomponent


def test_resp_napthalene_orca(tmp_path, monkeypatch):
    fixtures = Path(__file__).parent / "inputs" / "resp_classes"
    resp_folder = get_input_dir("resp_classes")
    # copy all stuffs in the resp_folder to tmp_path (current working directory)
    for item in os.listdir(resp_folder):
        s = os.path.join(resp_folder, item)
        d = os.path.join(tmp_path, item)
        if not os.path.isfile(s):
            continue 
        with open(s, 'rb') as fsrc, open(d, 'wb') as fdst:
            fdst.write(fsrc.read())
    startmulticomponent(["-f", "naphthalene_chcl3_resp.json"])

    try:
        import parmed as pmd 
    except ImportError:
        pytest.skip("parmed not installed; skipping test")

    # check the following:
    # 1. there is one and only one Cl- anion in this cube
    # 2. the atomic charge on "NAP" residue is 1.0
    # 3. the number of CL3 residue is close to 481
    assert os.path.exists("naphthalene-chcl3-resp.prmtop")
    structure = pmd.load_file("naphthalene-chcl3-resp.prmtop")

    cl_anions = [
        res
        for res in structure.residues
        if res.name.strip().upper() == "CL-"
    ]
    assert len(cl_anions) == 1

    nap_residue = next(
        (res for res in structure.residues if res.name.strip().upper() == "NAP"),
        None,
    )
    assert nap_residue is not None
    nap_charge = sum(atom.charge for atom in nap_residue.atoms)
    assert pytest.approx(1.0, rel=1e-3) == nap_charge

    cl3_residues = [
        res
        for res in structure.residues
        if res.name.strip().upper() == "CL3"
    ]
    assert np.allclose(len(cl3_residues), 481, rtol=0.05)
