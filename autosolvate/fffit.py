import os
import shutil
import numpy as np
import mdtraj as md

from .molecule import Molecule
from .utils import runMM, read_qm_folder
from .utils import PascalJobFileGenerator, SbatchJobManager
from .dockers import TleapDocker, MDGXDocker, ForceBalanceDocker


def run_mm(mol:Molecule, workingpath:str="./"):
    """
    Do conformational sampling using molecular mechanics
    """
    cwd = os.getcwd()
    if not os.path.exists(workingpath):
        os.makedirs(workingpath, exist_ok=True)
    shutil.copyfile(mol.prmtop, os.path.join(workingpath, os.path.basename(mol.prmtop)))
    shutil.copyfile(mol.inpcrd, os.path.join(workingpath, os.path.basename(mol.inpcrd)))
    name = os.path.basename(mol.prmtop)
    name = os.path.splitext(name)[0]

    for mminput in ["mmmin.in", "mmheat.in", "mmnvt.in", "mmnve.in"]:
        mminppath = os.path.join("templates", mminput)
        shutil.copyfile(mminppath, os.path.join(workingpath, mminput))
    os.chdir(workingpath)
    runMM(name)
    os.chdir(cwd)
    return

def run_qm(mol:Molecule, mmtrajectory:str, workingpath:str="./", nproc:int=16):
    """
    calculate QM energy and gradient
    """
    os.makedirs(workingpath, exist_ok=True)
    traj = md.load_netcdf(mmtrajectory, top=mol.prmtop)
    for i in range(traj.n_frames):
        t = traj.slice(i)
        t.save_xyz(os.path.join(workingpath, f"{mol.name}-{i}.xyz"))
    cwd = os.getcwd()
    os.chdir(workingpath)
    inputfiles = []
    inst = PascalJobFileGenerator()
    for i in range(traj.n_frames):
        jname = f"{mol.name}-{i}"
        inputfile = inst.gen_terachem_input_file(
            jobname=jname,
            run = "gradient",
            coord = jname+".xyz",
            method = "b3lyp",
            basis = "def2-svp",
            charge = 1,
            spinmult= 1,
            precision="dynamic",
            xtol = 1e-4,
        )
        inputfiles.append(inputfile)
    nrun_per_proc = int(np.ceil(len(inputfiles)/nproc))
    npprc = nrun_per_proc
    slurmscripts = []
    for i in range(nproc):
        scriptname = inst.gen_slurm_script_pascal(
            inputfile = inputfiles[i*npprc:(i+1)*npprc],
            jobname = f"fffit-{i}",
            excludenodelist=[1,2],
            software = "TeraChem/build2022",
        )
        slurmscripts.append(scriptname)
    jm = SbatchJobManager(slurmscripts, maxtime = 60*60)
    jm.wait_until_finish()
    os.chdir(cwd)


