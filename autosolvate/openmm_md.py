from __future__ import annotations

import os
from dataclasses import dataclass
from sys import stdout
from typing import Optional, Literal

from openmm import LangevinMiddleIntegrator
from openmm import XmlSerializer
from openmm import MonteCarloBarostat
from openmm import Platform
from openmm.app import AmberInpcrdFile, AmberPrmtopFile, Simulation
from openmm.app import HBonds
from openmm.app import DCDReporter, StateDataReporter, PDBReporter
from openmm.app import PME
from openmm.unit import kelvin, picosecond, femtoseconds, nanometer, bar


@dataclass(frozen=True)
class OpenMMRunOutputs:
    log_path: str
    dcd_path: str
    pdb_path: Optional[str]
    state_xml_path: str

def _require_file(path: str) -> str:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def _select_platform(prefer_gpu: bool = True) -> "object":
    if prefer_gpu:
        for name in ("CUDA", "OpenCL"):
            try:
                return Platform.getPlatformByName(name)
            except Exception:
                pass
    return Platform.getPlatformByName("CPU")

def append_reporters(
    simulation: Simulation,
    n_steps: int,
    report_interval: int = -1,
    snapshot_interval: int = -1,
    dcd_path: Optional[str] = None,
    log_path: Optional[str] = None,
    pdb_path: Optional[str] = None,
    stdout_reporter: bool = False,
    include_volume_density: bool = False,
):
    # preprocess parameters
    if report_interval <= 0:
        report_interval = min(100, max(1, n_steps // 100))
    if snapshot_interval <= 0:
        snapshot_interval = report_interval
    if dcd_path is not None:
        simulation.reporters.append(DCDReporter(dcd_path, report_interval))
    if log_path is not None:
        kwargs = dict(
            step=True,
            time=True,
            potentialEnergy=True,
            kineticEnergy=True,
            totalEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=n_steps,
            separator="\t",
        )
        if include_volume_density:
            kwargs["volume"] = True
            kwargs["density"] = True
        simulation.reporters.append(StateDataReporter(log_path, report_interval, **kwargs))
    if pdb_path is not None:
        simulation.reporters.append(PDBReporter(pdb_path, snapshot_interval))
    if stdout_reporter:
        stdout_kwargs = dict(
            step=True,
            time=True,
            potentialEnergy=True,
            temperature=True,
            progress=True,
            remainingTime=True,
            speed=True,
            totalSteps=n_steps,
            separator="\t",
        )
        if include_volume_density:
            stdout_kwargs["volume"] = True
            stdout_kwargs["density"] = True
        simulation.reporters.append(StateDataReporter(stdout, report_interval, **stdout_kwargs))



def run_openmm_minimize(
    *,
    prmtop_path: str,
    inpcrd_path: str,
    out_dir: str,
    temperature_k: float = 300.0,
    friction_per_ps: float = 1.0,
    timestep_fs: float = 2.0,
    nonbonded_cutoff_nm: float = 1.0,
    constraint_tolerance: float = 1e-5,
    platform_prefer_gpu: bool = True,
    max_iterations: int = 500,
) -> OpenMMRunOutputs:
    """Minimize an Amber-formatted system using OpenMM.

    Inputs must be Amber: `*.prmtop` and `*.inpcrd`.
    Writes a PDB snapshot, a serialized state XML, and a log.
    """

    os.makedirs(out_dir, exist_ok=True)

    prmtop = AmberPrmtopFile(prmtop_path)
    inpcrd = AmberInpcrdFile(inpcrd_path)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff_nm * nanometer,
        constraints=HBonds,
    )

    integrator = LangevinMiddleIntegrator(
        temperature_k * kelvin,
        friction_per_ps / picosecond,
        timestep_fs * femtoseconds,
    )
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = _select_platform(prefer_gpu=platform_prefer_gpu)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    log_path = os.path.join(out_dir, "openmm_minimize.log")
    dcd_path = os.path.join(out_dir, "openmm_minimize.dcd")
    pdb_path = os.path.join(out_dir, "openmm_minimize.pdb")
    state_xml_path = os.path.join(out_dir, "openmm_minimize_state.xml")

    with open(log_path, "w", encoding="utf-8") as f:
        f.write(f"platform={platform.getName()}\n")
        f.write(f"temperature_k={temperature_k}\n")
        f.write(f"timestep_fs={timestep_fs}\n")
        f.write("stage=minimize\n")
    simulation.minimizeEnergy(maxIterations=max_iterations)

    # Snapshot and state
    from openmm.app import PDBFile

    state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(pdb_path, "w", encoding="utf-8") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

    from openmm import XmlSerializer

    with open(state_xml_path, "w", encoding="utf-8") as f:
        f.write(XmlSerializer.serialize(state))

    return OpenMMRunOutputs(
        log_path=log_path,
        dcd_path=dcd_path,
        pdb_path=pdb_path,
        state_xml_path=state_xml_path,
    )


def run_openmm_nvt(
    *,
    prmtop_path: str,
    inpcrd_path: str,
    out_dir: str,
    n_steps: int = 1000,
    report_interval: int = 100,
    snapshot_interval: int = 500,
    temperature_k: float = 300.0,
    friction_per_ps: float = 1.0,
    timestep_fs: float = 2.0,
    nonbonded_cutoff_nm: float = 1.0,
    constraint_tolerance: float = 1e-5,
    platform_prefer_gpu: bool = True,
    starting_state_xml: Optional[str] = None,
    stdout_reporter: bool = True,
    write_pdb: bool = True,
    minimize_before_run: bool = True,
    minimize_max_iterations: int = 200,
) -> OpenMMRunOutputs:
    """Run an NVT simulation (Langevin).

    Outputs:
    - DCD trajectory (`openmm_prod.dcd`)
    - periodic PDB snapshots (`openmm_prod_stepXXXXX.pdb`)
    - a human-readable `StateDataReporter` log (`openmm_prod.log`)
    - final state XML (`openmm_prod_state.xml`)
    """

    os.makedirs(out_dir, exist_ok=True)

    prmtop = AmberPrmtopFile(prmtop_path)
    inpcrd = AmberInpcrdFile(inpcrd_path)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff_nm * nanometer,
        constraints=HBonds,
    )

    integrator = LangevinMiddleIntegrator(
        temperature_k * kelvin,
        friction_per_ps / picosecond,
        timestep_fs * femtoseconds,
    )
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = _select_platform(prefer_gpu=platform_prefer_gpu)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    if starting_state_xml is not None:
        with open(starting_state_xml, "r", encoding="utf-8") as f:
            state = XmlSerializer.deserialize(f.read())
        simulation.context.setState(state)

    # Ensure velocities are initialized; inpcrd often has positions only.
    simulation.context.setVelocitiesToTemperature(temperature_k * kelvin)

    log_path = os.path.join(out_dir, "openmm_nvt.log")
    dcd_path = os.path.join(out_dir, "openmm_nvt.dcd")
    pdb_path = os.path.join(out_dir, "openmm_nvt.pdb") if write_pdb else None
    state_xml_path = os.path.join(out_dir, "openmm_nvt_state.xml")

    append_reporters(
        simulation,
        n_steps=n_steps,
        report_interval=report_interval,
        snapshot_interval=snapshot_interval,
        dcd_path=dcd_path,
        log_path=log_path,
        pdb_path=pdb_path,
        stdout_reporter=stdout_reporter,
        include_volume_density=False,
    )

    if minimize_before_run:
        simulation.minimizeEnergy(maxIterations=minimize_max_iterations)

    # run with periodic snapshots
    simulation.step(n_steps)

    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(state_xml_path, "w", encoding="utf-8") as f:
        f.write(XmlSerializer.serialize(final_state))

    return OpenMMRunOutputs(
        log_path=log_path,
        dcd_path=dcd_path,
        pdb_path=pdb_path,
        state_xml_path=state_xml_path,
    )


def run_openmm_heat_nvt(
    *,
    prmtop_path: str,
    inpcrd_path: str,
    out_dir: str,
    n_steps: int = 1000,
    t_initial_k: float = 0.0,
    t_final_k: float = 300.0,
    friction_per_ps: float = 1.0,
    timestep_fs: float = 1.0,
    report_interval: int = 100,
    snapshot_interval: int = 500,
    nonbonded_cutoff_nm: float = 1.0,
    constraint_tolerance: float = 1e-5,
    platform_prefer_gpu: bool = True,
    starting_state_xml: Optional[str] = None,
    stdout_reporter: bool = True,
    write_pdb: bool = True,
    minimize_before_run: bool = True,
    minimize_max_iterations: int = 200,
) -> OpenMMRunOutputs:
    """Heat under NVT by ramping the bath temperature from `t_initial_k` to `t_final_k`.

    This mimics the classic "heat from 0 K" stage. We implement it by stepping in
    chunks and updating the integrator temperature linearly.
    """

    os.makedirs(out_dir, exist_ok=True)

    prmtop = AmberPrmtopFile(prmtop_path)
    inpcrd = AmberInpcrdFile(inpcrd_path)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff_nm * nanometer,
        constraints=HBonds,
    )

    integrator = LangevinMiddleIntegrator(
        max(t_initial_k, 0.0) * kelvin,
        friction_per_ps / picosecond,
        timestep_fs * femtoseconds,
    )
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = _select_platform(prefer_gpu=platform_prefer_gpu)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    if starting_state_xml is not None:
        with open(starting_state_xml, "r", encoding="utf-8") as f:
            state = XmlSerializer.deserialize(f.read())
        simulation.context.setState(state)

    # Initialize velocities at the initial temperature.
    simulation.context.setVelocitiesToTemperature(max(t_initial_k, 0.0) * kelvin)

    log_path = os.path.join(out_dir, "openmm_heat.log")
    dcd_path = os.path.join(out_dir, "openmm_heat.dcd")
    pdb_path = os.path.join(out_dir, "openmm_heat.pdb") if write_pdb else None
    state_xml_path = os.path.join(out_dir, "openmm_heat_state.xml")

    append_reporters(
        simulation,
        n_steps=n_steps,
        report_interval=report_interval,
        snapshot_interval=snapshot_interval,
        dcd_path=dcd_path,
        log_path=log_path,
        pdb_path=pdb_path,
        stdout_reporter=stdout_reporter,
        include_volume_density=False,
    ) 
    
    if minimize_before_run:
        simulation.minimizeEnergy(maxIterations=minimize_max_iterations)

    # Heat with linear ramp.
    steps_done = 0
    update_every = max(1, min(snapshot_interval, report_interval))
    while steps_done < n_steps:
        chunk = min(update_every, n_steps - steps_done)
        simulation.step(chunk)
        steps_done += chunk

        frac = steps_done / float(n_steps)
        target = t_initial_k + frac * (t_final_k - t_initial_k)
        integrator.setTemperature(target * kelvin)

    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(state_xml_path, "w", encoding="utf-8") as f:
        f.write(XmlSerializer.serialize(final_state))

    return OpenMMRunOutputs(
        log_path=log_path,
        dcd_path=dcd_path,
        pdb_path=pdb_path,
        state_xml_path=state_xml_path,
    )


def run_openmm_npt(
    *,
    prmtop_path: str,
    inpcrd_path: str,
    out_dir: str,
    n_steps: int = 1000,
    temperature_k: float = 300.0,
    pressure_bar: float = 1.0,
    barostat_interval: int = 25,
    friction_per_ps: float = 1.0,
    timestep_fs: float = 2.0,
    report_interval: int = 100,
    snapshot_interval: int = 500,
    nonbonded_cutoff_nm: float = 1.0,
    constraint_tolerance: float = 1e-5,
    platform_prefer_gpu: bool = True,
    starting_state_xml: Optional[str] = None,
    stdout_reporter: bool = True,
    write_pdb: bool = True,
    minimize_before_run: bool = True,
    minimize_max_iterations: int = 200,
) -> OpenMMRunOutputs:
    """Run an NPT simulation (Langevin + MonteCarloBarostat)."""

    os.makedirs(out_dir, exist_ok=True)

    prmtop = AmberPrmtopFile(prmtop_path)
    inpcrd = AmberInpcrdFile(inpcrd_path)

    system = prmtop.createSystem(
        nonbondedMethod=PME,
        nonbondedCutoff=nonbonded_cutoff_nm * nanometer,
        constraints=HBonds,
    )
    system.addForce(MonteCarloBarostat(pressure_bar * bar, temperature_k * kelvin, barostat_interval))

    integrator = LangevinMiddleIntegrator(
        temperature_k * kelvin,
        friction_per_ps / picosecond,
        timestep_fs * femtoseconds,
    )
    integrator.setConstraintTolerance(constraint_tolerance)

    platform = _select_platform(prefer_gpu=platform_prefer_gpu)
    simulation = Simulation(prmtop.topology, system, integrator, platform)
    simulation.context.setPositions(inpcrd.positions)
    if inpcrd.boxVectors is not None:
        simulation.context.setPeriodicBoxVectors(*inpcrd.boxVectors)

    if starting_state_xml is not None:
        with open(starting_state_xml, "r", encoding="utf-8") as f:
            state = XmlSerializer.deserialize(f.read())
        simulation.context.setState(state)

    simulation.context.setVelocitiesToTemperature(temperature_k * kelvin)

    log_path = os.path.join(out_dir, "openmm_npt.log")
    dcd_path = os.path.join(out_dir, "openmm_npt.dcd")
    pdb_path = os.path.join(out_dir, "openmm_npt.pdb") if write_pdb else None
    state_xml_path = os.path.join(out_dir, "openmm_npt_state.xml")

    append_reporters(
        simulation,
        n_steps=n_steps,
        report_interval=report_interval,
        snapshot_interval=snapshot_interval,
        dcd_path=dcd_path,
        log_path=log_path,
        pdb_path=pdb_path,
        stdout_reporter=stdout_reporter,
        include_volume_density=True,
    )

    if minimize_before_run:
        simulation.minimizeEnergy(maxIterations=minimize_max_iterations)

    simulation.step(n_steps)

    final_state = simulation.context.getState(getPositions=True, getEnergy=True)
    with open(state_xml_path, "w", encoding="utf-8") as f:
        f.write(XmlSerializer.serialize(final_state))

    return OpenMMRunOutputs(
        log_path=log_path,
        dcd_path=dcd_path,
        pdb_path=pdb_path,
        state_xml_path=state_xml_path,
    )
