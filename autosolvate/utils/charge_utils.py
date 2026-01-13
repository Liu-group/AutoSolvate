"""Utilities for inferring molecular charges from structure files.

This module mirrors the logic used in ``autoMCPB.py`` to derive a molecule's
formal charge from its structure. It is intended for closed-shell, main-group
molecules; transition metals and other atypical elements are not handled.
"""
import logging
import os
import subprocess
import tempfile
from typing import List, Optional

from rdkit import Chem
from rdkit.Chem import rdDetermineBonds
from rdkit import rdBase
from rdkit import RDLogger

from ..Common import OBABEL

# Allowed elements for automatic charge inference (upper-case).
_ALLOWED_MAIN_GROUP = {
    "H",
    "B",
    "C",
    "N",
    "O",
    "F",
    "P",
    "S",
    "CL",
    "BR",
    "I",
    "SI",
    "AL",
    "NA",
    "K",
    "LI",
    "MG",
    "CA",
}

_VALENCE_ELECTRONS = {
    "H": 1,
    "B": 3,
    "C": 4,
    "N": 5,
    "O": 6,
    "F": 7,
    "P": 5,
    "S": 6,
    "CL": 7,
    "BR": 7,
    "I": 7,
}

_VALENCE_SHELL = {
    "H": 2,
    "B": 8,
    "C": 8,
    "N": 8,
    "O": 8,
    "F": 8,
    "P": 8,
    "S": 8,
    "CL": 8,
    "BR": 8,
    "I": 8,
}


def _read_xyz_elements(xyz_path: str) -> List[str]:
    elements: List[str] = []
    with open(xyz_path, "r") as f:
        lines = f.readlines()
    for line in lines[2:]:
        parts = line.split()
        if len(parts) < 4:
            continue
        elements.append(parts[0].upper())
    return elements


def _convert_xyz_to_sdf(xyz_path: str, logger: logging.Logger) -> Optional[str]:
    """Convert xyz to a temporary sdf and return its path."""
    tmp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".sdf")
    tmp_file.close()
    sdf_path = tmp_file.name
    cmd = f"{OBABEL} -ixyz {xyz_path} -osdf -O {sdf_path} --errorlevel 0"
    try:
        subprocess.run(
            cmd,
            shell=True,
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
    except subprocess.CalledProcessError as exc:
        logger.warning("Auto charge: failed to convert %s to sdf (%s).", xyz_path, exc)
        # try:
        #     os.remove(sdf_path)
        # except OSError:
        #     pass
        return None
    return sdf_path


def _charge_by_enumeration(xyz_path: str, logger: logging.Logger) -> Optional[int]:
    """Enumerate possible charges and pick the one with minimal |atom formal charge| sum."""
    best_charge: Optional[int] = None
    best_abs_sum: Optional[float] = None
    for assumed in range(-5, 6):
        try:
            mol = Chem.MolFromXYZFile(xyz_path)
            if mol is None:
                continue
            rdDetermineBonds.DetermineConnectivity(mol)
            rdDetermineBonds.DetermineBondOrders(mol, charge=assumed)
            formal_charge = Chem.GetFormalCharge(mol)
            abs_sum = sum(abs(atom.GetFormalCharge()) for atom in mol.GetAtoms())
            if best_abs_sum is None or abs_sum < best_abs_sum:
                best_abs_sum = abs_sum
                best_charge = int(formal_charge)
        except Exception as exc:  # pragma: no cover - RDKit heuristics can fail
            logger.debug("Auto charge: failed at assumed charge %s (%s)", assumed, exc)
            continue
    if best_charge is not None:
        logger.info(
            "Auto charge: enumerated charges - selected %s with minimal |formal charges| sum %s",
            best_charge,
            best_abs_sum,
        )
    return best_charge


def _formal_charge_from_rdkit(sdf_path: str, logger: logging.Logger) -> Optional[int]:
    try:
        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
        mols = [mol for mol in suppl if mol is not None]
    except Exception as exc:  # pragma: no cover - extremely rare RDKit failure
        logger.warning("Auto charge: RDKit failed on %s (%s).", sdf_path, exc)
        return None
    if not mols:
        logger.warning("Auto charge: empty molecule read from %s.", sdf_path)
        return None
    mol = mols[0]
    try:
        total_charge = Chem.GetFormalCharge(mol)
    except Exception as exc:  # pragma: no cover - extremely rare RDKit failure
        logger.warning("Auto charge: RDKit cannot compute formal charge (%s).", exc)
        logger.warning("Auto charge: trying to sum individual atom formal charges.")
        try:
            total_charge = sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
        except Exception as exc:  # pragma: no cover - defensive
            logger.warning("Auto charge: cannot evaluate formal charge (%s).", exc)
            return None
    return int(total_charge)


def _valence_charge_from_sdf(sdf_path: str, logger: logging.Logger) -> int:
    """Fallback charge estimate using valence counting (mirrors autoMCPB)."""
    with open(sdf_path, "r") as f:
        data = f.readlines()
    atom_count = int(data[3][:3])
    bond_count = int(data[3][3:6])
    atomtypes = data[4 : 4 + atom_count]
    connectioninfo = data[4 + atom_count : bond_count + 4 + atom_count]

    linker_dic = {i: 0 for i in range(atom_count)}
    type_dic = {}
    connected: List[set] = []

    for idx, line in enumerate(atomtypes):
        if len(line.split()) > 3:
            type_dic[idx] = line.split()[3].upper()
    for line in connectioninfo:
        atomx = int(line[:3]) - 1
        atomy = int(line[3:6]) - 1
        bondorder = int(line[8])
        connected.append({atomx, atomy})
        linker_dic[atomx] += bondorder
        linker_dic[atomy] += bondorder

    ligand_charge = 0
    for atom, atype in type_dic.items():
        if atype not in _VALENCE_ELECTRONS:
            logger.warning("Auto charge: element %s not supported in valence fallback.", atype)
            return 0
        valence_electron = _VALENCE_ELECTRONS[atype]
        longpair = _VALENCE_SHELL[atype] - (int(valence_electron) + linker_dic[atom])
        if longpair < 0:
            ligand_charge += abs(longpair)
        else:
            ligand_charge -= abs(longpair)
    return int(ligand_charge)


def infer_charge_from_xyz(xyz_path: str, logger: Optional[logging.Logger] = None) -> int:
    """
    Infer the integer charge for a main-group molecule from an xyz file.

    The routine first checks for unsupported elements, then relies on RDKit
    formal charges (after obabel conversion). If RDKit cannot parse the
    structure, it falls back to a valence-electron counting heuristic similar
    to the one in ``autoMCPB.py``.

    Parameters
    ----------
    xyz_path: str
        Path to the xyz file describing the molecule.
    logger: logging.Logger, optional
        Logger to report warnings/info; a module logger is used if omitted.

    Returns
    -------
    int
        The inferred integer charge. Returns 0 when inference is not possible.
    """
    RDLogger.DisableLog('*')
    log = logger if logger is not None else logging.getLogger(__name__)

    elements = _read_xyz_elements(xyz_path)
    for elem in elements:
        if elem.upper() not in _ALLOWED_MAIN_GROUP:
            log.warning(
                "Auto charge: element %s is not supported; default charge 0 will be used.",
                elem,
            )
            return 0

    # First try enumeration within [-5, 5] using bond-order determination.
    enum_charge = _charge_by_enumeration(xyz_path, log)
    if enum_charge is not None:
        return int(enum_charge)

    sdf_path = _convert_xyz_to_sdf(xyz_path, log)
    if sdf_path is None:
        return 0

    charge = _formal_charge_from_rdkit(sdf_path, log)
    if charge is None:
        charge = _valence_charge_from_sdf(sdf_path, log)

    try:
        os.remove(sdf_path)
    except OSError:
        pass
    RDLogger.EnableLog('*')
    return int(charge)
