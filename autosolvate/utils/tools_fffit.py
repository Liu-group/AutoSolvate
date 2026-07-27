#--------------------------------------------------------------------------------------------------#
# Utilaity functions for fitting force field parameters
# author: Fangning Ren (2023-04-27)
# path: autosolvate/utils/tools_fffit.py
#--------------------------------------------------------------------------------------------------#

import os
import re
import subprocess
from typing import Iterable, List, Tuple, Union, TextIO
import numpy as np

from .tools import read_xyz

def read_gradient(fname:str):
    gradientdata = []
    start = False
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("Gradient units are Hartree/Bohr"):
                start = True
            if line.startswith("Net gradient:"):
                start = False
            if start == True:
                gradientdata.append(line)

    gradientdata = gradientdata[3:-1]
    natom = len(gradientdata)
    gradient = np.zeros((natom, 3), dtype = float)
    for i in range(natom):
        # print(gradientdata[i])
        gradient[i] = list(map(float, gradientdata[i].split()))
    return gradient

def read_energy(fname:str):
    eline = "FINAL ENERGY: 0.000000000 a.u."
    with open(fname, "r") as f:
        for line in f:
            if line.startswith("FINAL ENERGY:"):
                eline = line
    eline = eline.replace("FINAL ENERGY:", "").replace("a.u.", "")
    return float(eline)

def read_qm_folder(folder:str, prefix:str):
    fnames = os.listdir(folder)
    fnames = [fname for fname in fnames if fname.startswith(prefix)]
    cnames = [os.path.splitext(fname)[0] for fname in fnames]
    cnames = sorted(list(set(cnames)))
    xyzs, grds, enes = [], [], []
    for cname in cnames:
        xyzpath = os.path.join(folder, cname+".xyz")
        outpath = os.path.join(folder, cname+".out")
        if not os.path.exists(xyzpath) or not os.path.exists(outpath):
            continue
        try:
            elem, coord = read_xyz(xyzpath)
            grad = read_gradient(outpath)
            ener = read_energy(outpath)
            if coord.shape != grad.shape:
                raise Exception("coord and grad shape mismatch")
            if ener == 0.0:
                raise Exception("energy is zero")
            xyzs.append(coord)
            grds.append(grad)
            enes.append(ener)
        except Exception as e:
            print("Error when reading data from", xyzpath)
            print(e)
            continue
    return xyzs, grds, enes

def write_mdgx_dat(fname:str, energys):
    with open(fname, "w") as f:
        for e in energys:
            f.write("{:>16.8f}\n".format(e))

def get_mdgx_refdat(folder:str, prefix:str, target:str, nconf:int):
    folderpath = folder
    energys = []
    for i in range(nconf):
        outpath = f"{folderpath}/{prefix}-{i}.out"
        energy = read_energy(outpath)
        print(f"reading frame {i}. energy {energy}")
        energys.append(energy)
    write_mdgx_dat(target, prefix+".netcdf", prefix+".prmtop", energys)
    return energys

def runMM(filename='water_solvated', freeze_solute=False):        
   print('MM Energy minimization')
   cmd=' -O -i mmmin.in -o mmmin.out -p '+filename+'.prmtop -c '+filename+'.inpcrd -r mm.ncrst -inf mmmin.info'
   if freeze_solute:
      cmd = cmd + ' -ref '+filename+'.inpcrd '
   cmd= 'sander'+ cmd
   subprocess.call(cmd, shell=True)
   print('MM Heating')
   cmd=' -O -i mmheat.in -o mmheat.out -p '+filename+'.prmtop -c mm.ncrst -r mm.ncrst -x '+filename+'-heat.netcdf -inf mmheat.info'
   if freeze_solute:
      cmd = cmd + ' -ref '+filename+'.inpcrd '
   cmd= 'sander'+ cmd
   subprocess.call(cmd, shell=True)
   print('MM NVE equilibration')
   cmd=' -O -i mmnve.in -o mmnve.out -p '+filename+'.prmtop -c mm.ncrst -r mm.ncrst -x '+filename+'-mmnve.netcdf -inf mmnve.info'
   if freeze_solute:
      cmd = cmd + ' -ref '+filename+'.inpcrd '
   cmd= 'sander'+ cmd
   subprocess.call(cmd, shell=True)
   print('MM NVT equilibration')
   cmd=' -O -i mmnvt.in -o mmnvt.out -p '+filename+'.prmtop -c mm.ncrst -r mm.ncrst -x '+filename+'-mmnvt.netcdf -inf mmnvt.info'
   if freeze_solute:
      cmd = cmd + ' -ref '+filename+'.inpcrd '
   cmd= 'sander'+ cmd
   subprocess.call(cmd, shell=True)

class PascalJobFileGenerator:
    # This class only works for our cluster
    # This class also has duplicate functionality with TeraChemDocker
    # TODO: make it more general
    def __init__(self, cluster = "pascal"):
        if cluster == "pascal":
            self.gen_slurm_script = self.gen_slurm_script_pascal
        elif cluster == "bridges2":
            self.gen_slurm_script = self.gen_slurm_script_bridges2

    def gen_slurm_script_bridges2(self, inputfile,
        jobname = None,
        software = "TeraChem/mpich2",
        partition = "day-long", 
        nproc = 1, 
        gpu = True,
        excludenodelist = []
        ):
        if isinstance(inputfile, str):
            jobname = inputfile.replace(".inp", "")
        elif isinstance(inputfile, Iterable):
            jobname = "_".join(inputfile) if jobname == None else jobname
        subfile = jobname + ".sh"
        myfile = open(subfile,'w')
        myfile.write("#!/bin/bash\n")
        if partition == "short":
            myfile.write("#SBATCH -t 00:30:00\n")
        elif partition == "day-long":
            myfile.write("#SBATCH -t 24:00:00\n")
        elif partition == "week-long":
            myfile.write("#SBATCH -t 7-00:00:00\n")
        elif partition == "month-long":
            myfile.write("#SBATCH -t 30-00:00:00\n")
        
        myfile.write("#SBATCH -N 1\n")
        myfile.write("#SBATCH -p GPU-shared\n")
        myfile.write("#SBATCH --gpus=1\n")
        myfile.write("#SBATCH -A che210029p")
        myfile.write(f"#SBATCH --job-name={jobname}")
        myfile.write("set -x\n")
        myfile.write("echo $LD_LIBRARY_PATH\n")
        myfile.write("module use /jet/home/ffangliu/modulefiles\n")
        myfile.write("module load terachem/v02082021\n")
        
        if isinstance(inputfile, str):
            myfile.write(f"terachem {inputfile} > {jobname}.out\n")
            myfile.write(f"cp scr-{jobname}/{jobname}.molden ./\n")
        elif isinstance(inputfile, Iterable):
            for ifile in inputfile:
                ifile = ifile.replace(".inp", "")
                myfile.write(f"terachem {ifile}.inp > {ifile}.out\n")
                myfile.write(f"cp scr-{ifile}/{ifile}.molden ./\n")

    def gen_slurm_script_pascal(self, inputfile,
        jobname = None,
        software = "TeraChem/mpich2",
        partition = "day-long", 
        nproc = 1, 
        gpu = True,
        copy_molden = False,
        copy_optimized_xyz = False,
        excludenodelist = []
        ):
        if isinstance(inputfile, str):
            jobname = inputfile.replace(".inp", "")
        elif isinstance(inputfile, Iterable):
            jobname = "_".join(inputfile) if jobname == None else jobname
        subfile = jobname + ".sh"
        myfile = open(subfile,'w')
        myfile.write("#!/bin/bash\n")
        if partition == "short":
            myfile.write("#SBATCH --time=00:30:00\n")
            myfile.write("#SBATCH --partition=short\n")
        elif partition == "day-long":
            myfile.write("#SBATCH --time=24:00:00\n")
            myfile.write("#SBATCH --partition=day-long\n")
        elif partition == "week-long":
            myfile.write("#SBATCH --time=7-00:00:00\n")
            myfile.write("#SBATCH --partition=week-long\n")
        elif partition == "month-long":
            myfile.write("#SBATCH --time=30-00:00:00\n")
            myfile.write("#SBATCH --partition=month-long\n")
        if nproc == 1:
            myfile.write("#SBATCH --nodes=1\n")
            myfile.write("#SBATCH --mem=1G\n")
            myfile.write("#SBATCH --ntasks=1\n")
            myfile.write("#SBATCH --cpus-per-task=1\n")
        else:
            myfile.write("#SBATCH --nodes=1\n")
            myfile.write("#SBATCH --mem=%dG\n" %(nproc * 4))
            myfile.write("#SBATCH --ntasks=%d\n" %(nproc))
            myfile.write("#SBATCH --cpus-per-task=1\n")
        if gpu:
            myfile.write("#SBATCH --gres=gpu:1\n")
        if excludenodelist:
            myfile.write("#SBATCH --exclude=node[" + ",".join([str(i) for i in excludenodelist]) + "]\n")

        if software.lower().find("terachem") != -1:
            myfile.write("echo $HOSTNAME\n")
            myfile.write("module load %s\n" %software)
            myfile.write("echo $CUDA_VISIBLE_DEVICES\n\n")
            if isinstance(inputfile, str):
                inputfile = [inputfile,]
            if isinstance(inputfile, Iterable):
                for ifile in inputfile:
                    ifile = ifile.replace(".inp", "")
                    myfile.write(f"terachem {ifile}.inp > {ifile}.out\n")
                    if copy_molden:
                        myfile.write(f"cp scr-{ifile}/{ifile}.molden ./\n")
                    if copy_optimized_xyz:
                        newxyzname = ifile + "-optim.xyz"
                        myfile.write(f"mv scr-{ifile}/optim.xyz ./{newxyzname}\n")
            
        elif software.lower().find("orca") != -1:
            myfile.write(f"module load {software}\n")  
            myfile.write("cd $SLURM_SUBMIT_DIR\n")                                                        
            myfile.write('export PATH="//usr/local/mpi/gcc/openmpi-cuda11.3-4.1.1/bin:$PATH"\n')                                      
            if isinstance(inputfile, str):
                myfile.write(f"//opt/orca/5.0.2/orca {inputfile} > {jobname}.out\n")
                myfile.write(f"orca_2mkl -molden {jobname}")
            elif isinstance(inputfile, Iterable):
                for ifile in inputfile:
                    ifile = ifile.replace(".inp", "")
                    myfile.write(f"//opt/orca/5.0.2/orca {ifile}.inp > {ifile}.out\n")
                    myfile.write(f"orca_2mkl -molden {ifile}")
            
            myfile.write(f"//opt/orca/5.0.2/orca {inputfile} > {jobname}.out\n")
        myfile.close()
        return subfile

    def gen_terachem_input_file(self, jobname, coord,
        run = "energy",
        method = "rhf",
        basis = "6-31g",
        charge = "0",
        spinmult = "1",
        guess = "generate",
        maxit = "200",
        scf = "diis",
        precision = "double",
        levelshift = "no",
        levelshiftvala = "0.1",
        levelshiftvalb = "0.1",
        printitermin = "1",
        printitermax = "1",
        dispersion = "d3",
        rc_w = "",
        cis = "no",
        tddft = "no",
        cisnumstates = "2",
        cistarget = "1",
        cistransdensity = "no",
        cisdiffdensity = "no",
        epsilon = "1.0",
        fast_epsilon = "1.0",
        lr_pcm_solvation = "neq",
        xtol = "1.0e-6",
        convthre = "3.0e-5",
        sphericalbasis = "yes",
        additionalinfo = ""
    ):
        fnames = os.listdir(os.getcwd())
        if not os.path.exists(coord):
            coord = coord.replace(".xyz", ".pdb")
            if not os.path.exists("coord"):
                raise FileNotFoundError(f"cannot find {coord} or {coord.replace('.pdb', '.xyz')}.")
        jobname = jobname.replace(".inp", "")
        myfile = open(jobname+".inp", "w")
        # necessary information
        myfile.write("jobname        " + jobname + "\n")
        myfile.write("coordinates    " + coord + "\n")
        myfile.write("run            " + run + "\n")
        myfile.write("scrdir         scr-" + jobname + "\n")
        myfile.write("gpus           1\n")

        # necessary information for calculation
        if method in ["gfnxtb","gfn2xtb"]:  # ust semi emperical method
            basis = method
            guess = "hcore"
        myfile.write("method         " + str(method) + "\n")
        myfile.write("basis          " + str(basis) + "\n") 
        myfile.write("charge         " + str(charge) + "\n")
        myfile.write("spinmult       " + str(spinmult) + "\n")
        myfile.write("sphericalbasis " + str(sphericalbasis) + "\n")

        # scf setting
        myfile.write("guess          " + str(guess) + "\n")
        myfile.write("maxit          " + str(maxit) + "\n")
        myfile.write("scf            " + str(scf) + "\n")
        myfile.write("precision      " + str(precision) + "\n")
        myfile.write("convthre       " + str(convthre) + "\n")
        myfile.write("xtol           " + str(xtol) + "\n")
        if levelshift == "yes":
            myfile.write("levelshift     " + str(levelshift) + "\n")
            myfile.write("levelshiftvala " + str(levelshiftvala) + "\n")
            myfile.write("levelshiftvalb " + str(levelshiftvalb) + "\n")
        if int(printitermax) > int(printitermin):
            myfile.write("printitermin   " + str(printitermin) + "\n")
            myfile.write("printitermax   " + str(printitermax) + "\n")

        # dft setting
        myfile.write("dftgrid        1\n")
        myfile.write("dynamicgrid    no\n")
        if method.find("b97") != -1 and method != "wb97xd3":
            dispersion = "no"
        myfile.write("dispersion     " + str(dispersion) + "\n")
        if "w" in method and rc_w != "":    # omega tuning
            myfile.write("rc_w           " + str(rc_w) + "\n")

        # tddft setting
        if cis == "yes" or cis == True or tddft == "yes" or tddft == True:
            myfile.write("cis            yes\n")
            myfile.write("cisnumstates   " + str(cisnumstates) + "\n")
            myfile.write("cismult        " + str(spinmult) + "\n")
            myfile.write("cistarget      " + str(cistarget) + "\n")
            myfile.write("cismaxiter     100\n")
            myfile.write("cismax         " + str(int(cisnumstates) * 5) + "\n")
            myfile.write("cisprintthresh 0.01\n")
            myfile.write("cistransdensity        " + str(cistransdensity) + "\n")
            myfile.write("cisdiffdensity         " + str(cisdiffdensity) + "\n")

        # geometry optimization setting
        if run == "minimize":
            myfile.write("min_coordinates   tric\n")
            myfile.write("new_minimizer  yes\n")
            myfile.write("min_converge_e    5.0e-7\n")  # release the convergence threshold. other settings remains the same.
            myfile.write("purify         no\n")
            myfile.write("nstep          10000\n")

        # pcm setting
        dielec = float(str(epsilon))
        fielec = float(str(fast_epsilon))
        if dielec > 1.0:
            myfile.write("\npcm             cosmo\n")
            myfile.write("pcm_grid        iswig\n")
            myfile.write("pcmgrid_heavy   17\n")
            myfile.write("pcmgrid_h       17\n")
            myfile.write("epsilon         " + str(epsilon) + "\n")
            myfile.write("pcm_radii       bondi\n")
            myfile.write("dynamiccg       0\n")
        if dielec > 1.0 and fielec > 1.0 and (cis == "yes" or cis == True or tddft == "yes" or tddft == True):
            myfile.write("lr_pcm_solvation      " + str(lr_pcm_solvation) + "\n")
            myfile.write("fast_epsilon    " + str(fast_epsilon) + "\n")
        if dielec > 1.0 and fielec <=1.0 and (cis == "yes" or cis == True or tddft == "yes" or tddft == True):
            raise ValueError(f"fast epsilon {fielec} should be larger than 1.")

        myfile.write(additionalinfo)
        myfile.write("\nend\n")
        myfile.close()
        return jobname.replace(".inp", "")+".inp"

    def submit_job(self, submitfile):
        command = f"sbatch {submitfile}"
        subprocess.run(command, shell = True, cwd = os.getcwd())

    def run_terachem(self, inputfilename, build = "mpich2"):
        fname = inputfilename.replace(".inp", "")
        commands = [
            f"module load TeraChem/{build}",
            f"terachem {fname}.inp > {fname}.out",
            f"cp scr-{fname}/{fname}.molden ./",
            ]
        for command in commands:
            print(command)
            subprocess.run(command, shell = True)
