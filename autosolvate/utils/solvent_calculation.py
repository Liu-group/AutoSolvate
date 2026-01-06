"""Solvent and mixture calculation utilities.

Separated from tools.py to avoid circular imports and keep responsibilities clear.
"""
from __future__ import annotations

import math
from typing import Iterable, List, Sequence, Union

from ..Common import N_A, SOLVENT_DENSITY, SOLVENT_MW
from .tools import determine_mw_from_pdb, determine_mw_from_xyz

CubeSize = Union[float, Sequence[float]]


def calculate_solvent_number(solvent_name: str, volume_m3: float) -> int:
    """Calculate number of molecules for a pure solvent in a given volume (m^3)."""
    density = SOLVENT_DENSITY[solvent_name]
    weight = SOLVENT_MW[solvent_name]
    mass = volume_m3 * density  # kg
    mol = mass * 1000 / weight
    number = mol * N_A
    return int(number)


def calculate_solvent_number_from_density(molecular_weight_g_mol: float, density_g_cm3: float, volume_m3: float) -> int:
    """Calculate molecule count from density (g/cm^3), MW (g/mol), and volume (m^3)."""
    volume_cm3 = volume_m3 * 1e6
    mass_g = density_g_cm3 * volume_cm3
    mol = mass_g / molecular_weight_g_mol
    return int(mol * N_A)


def calculate_solvent_numbers_from_weight_portions(
    solvent_mws: Iterable[float],
    solvent_densities: Iterable[float],
    weight_portions: Iterable[float],
    total_volume_m3: float,
) -> List[int]:
    # Densities in input JSON are typically in g/cm^3. The mixture rule below expects
    # consistent units; we keep density in g/cm^3 and convert volume to cm^3.
    densities_g_cm3 = [float(d) if float(d) < 50 else float(d) / 1000.0 for d in solvent_densities]
    portions = [float(p) for p in weight_portions]
    mws_g_mol = [float(mw) for mw in solvent_mws]

    estimated_density_g_cm3 = 1.0 / sum(p / rho for p, rho in zip(portions, densities_g_cm3))
    volume_cm3 = float(total_volume_m3) * 1e6
    total_mass_g = volume_cm3 * estimated_density_g_cm3
    numbers = [total_mass_g * p / mw * N_A for mw, p in zip(mws_g_mol, portions)]
    return [int(n) for n in numbers]


def calculate_solvent_numbers_from_volume_portions(
    solvent_mws: Iterable[float],
    solvent_densities: Iterable[float],
    volume_portions: Iterable[float],
    total_volume_m3: float,
) -> List[int]:
    numbers = [portion * total_volume_m3 * density * 1000 * N_A * 1000 / mw for mw, density, portion in zip(solvent_mws, solvent_densities, volume_portions)]
    return list(map(int, numbers))


def calculate_solvent_numbers_from_molar_portions(
    solvent_mws: Iterable[float],
    solvent_densities: Iterable[float],
    molar_portions: Iterable[float],
    total_volume_m3: float,
) -> List[int]:
    weight_portions = [molar_portion * mw for molar_portion, mw in zip(molar_portions, solvent_mws)]
    total_weight_portion = sum(weight_portions)
    weight_portions = [wp / total_weight_portion for wp in weight_portions]
    estimated_density = 1 / sum(portion / density for portion, density in zip(weight_portions, solvent_densities))
    mass = total_volume_m3 * estimated_density  # kg
    mass_g = mass * 1000
    numbers = [mass_g * portion / mw * N_A for mw, portion in zip(solvent_mws, weight_portions)]
    return list(map(int, numbers))


def cube_size_to_volume_m3(cube_size: CubeSize) -> float:
    """Convert cube size (Angstrom, scalar or iterable of 3) to volume in m^3."""
    if isinstance(cube_size, (list, tuple)):
        if len(cube_size) != 3:
            raise ValueError("cube_size list/tuple must have 3 elements")
        a, b, c = map(float, cube_size)
        volume_ang3 = a * b * c
    else:
        volume_ang3 = float(cube_size) ** 3
    return volume_ang3 * 1e-30


def estimate_density_g_cm3(total_mass_g: float, volume_m3: float) -> float:
    """Return density in g/cm^3 given mass in g and volume in m^3."""
    volume_cm3 = volume_m3 * 1e6
    if volume_cm3 <= 0:
        raise ValueError("Volume must be positive")
    return total_mass_g / volume_cm3


def estimate_system_density(solutes: Sequence[dict], solvents: Sequence[dict], cube_size: CubeSize) -> float:
    """Estimate overall density (g/cm^3) from solute/solvent counts and MW."""
    volume_m3 = cube_size_to_volume_m3(cube_size)
    total_mass_g = 0.0
    for item in list(solutes) + list(solvents):
        if "number" not in item:
            number = 1
        else:
            number = int(item["number"])
        if "molecular_weight" not in item:
            if "xyzfile" in item:
                if item["xyzfile"].endswith(".pdb"):
                    mw = determine_mw_from_pdb(item["xyzfile"])
                else:
                    mw = determine_mw_from_xyz(item["xyzfile"])
            else:
                raise ValueError("Molecular weight must be provided if no xyzfile is given")
        else:
            mw = float(item["molecular_weight"])
        total_mass_g += number * mw / N_A
    return estimate_density_g_cm3(total_mass_g, volume_m3)


def estimate_solvent_density(solvents: Sequence[dict], cube_size: CubeSize) -> float:
    """Estimate density using only solvent components."""
    volume_m3 = cube_size_to_volume_m3(cube_size)
    total_mass_g = 0.0
    for item in solvents:
        if "number" not in item or "molecular_weight" not in item:
            continue
        try:
            number = int(item["number"])
            mw = float(item["molecular_weight"])
        except Exception:
            continue
        if number <= 0 or mw <= 0:
            continue
        total_mass_g += number * mw / N_A
    return estimate_density_g_cm3(total_mass_g, volume_m3)


def adjust_cube_size_for_density(current_cube_size: CubeSize, current_density: float, target_density: float) -> CubeSize:
    """Suggest a new cube size scaled to reach target density (g/cm^3)."""
    if current_density <= 0 or target_density <= 0:
        raise ValueError("Densities must be positive")
    scale = (current_density / target_density) ** (1.0 / 3.0)
    if isinstance(current_cube_size, (list, tuple)):
        return [float(x) * scale for x in current_cube_size]
    return float(current_cube_size) * scale


def scale_solvent_numbers(solvents: Sequence[dict], current_density: float, target_density: float) -> List[int]:
    """Scale solvent numbers uniformly to move density toward a target while keeping ratios."""
    if current_density <= 0:
        raise ValueError("Current density must be positive")
    factor = target_density / current_density
    # Use rounding to nearest int while keeping ratios approximately
    return [max(1, int(round(solvent.get("number", 0) * factor))) for solvent in solvents]
