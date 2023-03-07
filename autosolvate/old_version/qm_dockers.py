import sys
import os
import subprocess
import logging
import time
from abc import * 
from typing import List, Tuple, TextIO, Iterable, Union, Optional, Dict, Any

import numpy as np

from Common import *
from molecule import System, Molecule
from complex import MoleculeComplex
from solvated import SolvatedSystem
from dockers import GeneralDocker
from tools import *


class SbatchJob():
    def __init__(self, script:str, maxtime = 86400*7, submit_dir = "", check_interaval = 10, logger = None):
        """A data class that stores a slurm job information"""
        if not isinstance(script, str):
            raise ValueError("script file path required. the input is a "+ str(type(script)))
        if not os.path.exists(script):
            raise OSError("file not found: " + script)
        self.name = script
        with open(script, "r") as f:
            self.script = f.read()
        self.jobid          = -1
        self.partition      = ""
        self.time           = 0
        self.jobname        = ""
        self.user           = ""
        self.status         = "U"
        self.nodelist       = ""
        self.output         = ""
        if not submit_dir:
            submit_dir = os.path.dirname(script)
            if not submit_dir:
                submit_dir = os.getcwd()
        self.submit_dir = submit_dir

        self.maximum_time = maxtime
        self.check_interval = check_interaval

        if not isinstance(logger, logging.Logger):
            self.logger = logging.getLogger("SbatchJob")
            self.logger.setLevel(logging.INFO)
            self.logger.addHandler(logging.StreamHandler(sys.stdout))
        else:
            self.logger = logger

    def asctime2second(self, s:str):
        if s.count(":") == 1:
            sec, miu = int(s[-2:]), int(s[:-3])
            self.time = sec + 60 * miu
        elif s.count(":") == 2 and s.count("-") == 0:
            sec, miu, hor = int(s[-2:]), int(s[-5:-3]), int(s[:-6])
            self.time = sec + 60 * miu + 3600 * hor
        elif s.count(":") == 2 and s.count("-") == 1:
            d, s = s.split("-")
            sec, miu, hor = int(s[-2:]), int(s[-5:-3]), int(s[:-6])
            self.time = sec + 60 * miu + 3600 * hor + 86400 * int(d)
        else:
            self.logger.warn("warning: the time " + s + " cannot be processed.")
            self.time = self.time
        return self.time

    def check_finished(self):
        data = subprocess.getoutput("squeue").splitlines()
        if len(data) == 1:
            self.status = "F"
            return 
        for line in data[1:]:
            args = line.split()
            jid = int(args[0])
            if jid == self.jobid:
                self.partition  = args[1]
                self.jobname    = args[2]
                self.user       = args[3]
                self.status     = args[4]
                self.time       = self.asctime2second(args[5])
                self.nodelist   = args[-1]
                break
        else:
            self.status = "F"

    def submit(self):
        cwd = os.getcwd()
        if not os.path.exists(self.submit_dir):
            os.makedirs(self.submit_dir)
        os.chdir(self.submit_dir)
        result = subprocess.getoutput("sbatch " + self.name)
        os.chdir(cwd)
        jobid = result.split()[-1]
        try:
            self.jobid = int(jobid)
            self.output = os.path.join(os.getcwd(), "slurm-"+jobid+".out")
            self.check_finished()
        except Exception:
            self.logger.warn("failed to submit script " + self.name + ". slurm report:")
            self.logger.warn(result)

    def wait_until_finish(self):
        if self.status == "U":
            self.submit()
        for i in range(self.maximum_time // self.check_interval):
            time.sleep(self.check_interval)
            self.check_finished()
            if self.status == "F":
                self.logger.info(f"job {self.jobid:9d} {self.jobname:9s} finished. ")
                break
            else:
                self.logger.info(f"job {self.jobid:9d} {self.jobname:9s} status {self.status:2s} time {self.time:9d}")
        else:
            self.logger.warn(f"job {self.jobid:9d} {self.jobname:9s} has been running for more than {self.maximum_time} seconds, no longer waiting for this task to complete.")






class TeraChemDocker(GeneralDocker):
    """
    TeraChem Python Docker
    """
    _SUPPORT_INPUT_FORMAT = ["xyz", "pdb"]

    def __init__(self, 
                 terachem_path:         str = "",
                 jobname:               str = "",
                 workfolder:            str = WORKING_DIR,
                 exeoutfile:            str = None,
                 srun_use:              bool = False,
                 sbatch_use:            bool = True,
                 **kwargs
    ) -> None:
        super().__init__(
            executable = terachem_path,
            workfolder = workfolder,
            exeoutfile = exeoutfile)
        self.terachem       = self.get_tc_path(terachem_path)
        self.jobname        = jobname
        self.executable     = self.terachem
        self.terachem_out   = exeoutfile
        self.srun_use       = srun_use
        self.sbatch_use     = sbatch_use
        self.logger.name    = self.__class__.__name__

        self.coord_format   = "xyz"
        self.terachem_package_path = "TeraChem"

        self.input_file     = ""
        self.slurm_script   = ""
        self.output_file    = ""
        self.scr_path       = ""

        self.tc_args        = kwargs
        for key, value in kwargs.items():
            self.logger.info("Set {} to {}".format(key, value))
            self.tc_args[key] = value

    def get_tc_path(self, terachem_path:str) -> str:
        traversed_path = []
        traversed_path.append(terachem_path)
        if terachem_path == "" or not os.path.exists(terachem_path):
            terachem_path = "terachem"
            traversed_path.append(terachem_path)
            result = subprocess.getoutput("which terachem")
            self.terachem_package_path = "TeraChem"
            if result.find("no terachem in") != -1:
                tcmodules = subprocess.getoutput("module avail TeraChem")
                tcmodules = [s.replace("(default)", "", 1) for s in tcmodules.split() if s.find("TeraChem") != -1]
                res = subprocess.getoutput("module load TeraChem")
                self.terachem_package_path = "TeraChem"
                if res.find("ERROR") != -1:
                    for tcmodule in tcmodules:
                        res = subprocess.getoutput("module load {tcmodule}".format(tcmodule = tcmodules[0]))
                        self.terachem_package_path = tcmodules[0]
                        traversed_path.append(tcmodule)
                        if res.find("ERROR") == -1:
                            break
                    else:
                        self.logger.critical("No TeraChem module found in :{}".format(" ".join(traversed_path)))
                        raise ModuleNotFoundError("No terachem in: " + " ".join(traversed_path))
        self.terachem = terachem_path
        return terachem_path

    def check_system(self, mol:System):
        for fmd in TeraChemDocker._SUPPORT_INPUT_FORMAT:
            if mol.check_exist(fmd):
                self.logger.info("TeraChemDocker supports {fmd} format.".format(fmd = fmd)) 
                self.coord_format = fmd
                break
        else:
            self.logger.critical("TeraChemDocker does not support any of the following format: {fmd}".format(fmd = " ".join(TeraChemDocker._SUPPORT_INPUT_FORMAT)))
            raise ValueError("TeraChemDocker does not support any of the following format: {fmd}".format(fmd = " ".join(TeraChemDocker._SUPPORT_INPUT_FORMAT)))
        self.jobname = mol.name if self.jobname == "" else self.jobname
        self.logger.info("Checking input parameters...")
        charge, spinmult = 0, 1
        if isinstance(mol, Molecule):
            charge = mol.charge
        elif isinstance(mol, MoleculeComplex) or isinstance(mol, SolvatedSystem):
            charge = mol.netcharge
        spinmult = mol.multiplicity
        self.logger.info("\t jobname: {jobname}".format(jobname = self.jobname))
        self.logger.info("\t coord: {coord}".format(coord = getattr(mol, self.coord_format)))
        self.logger.info("\t charge: {charge}".format(charge = charge))
        self.logger.info("\t spinmult: {spinmult}".format(spinmult = spinmult))
        if "run" not in self.tc_args:
            self.tc_args["run"] = "energy"
        for key, value in self.tc_args.items():
            self.logger.info("\t {key}: {value}".format(key = key, value = value))
    
    def generate_tc_input(self, jobname:str, coord:str, charge:int, spinmult:int,
        run = "energy",
        method = "rhf",
        basis = "6-31g",
        guess = "generate",
        maxit = "200",
        scf = "diis",
        precision = "dynamic",
        levelshift = "no",
        levelshiftvala = "0.1",
        levelshiftvalb = "0.1",
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
        resp = "no",
        additionalinfo = ""
    ):
        fnames = os.listdir(os.getcwd())
        if not os.path.exists(coord):
            coord = coord.replace(".xyz", ".pdb")
            if not os.path.exists("coord"):
                raise FileNotFoundError(f"cannot find {coord} or {coord.replace('.pdb', '.xyz')}.")
        jobname = jobname.replace(".inp", "")
        inputfilename = jobname.replace(".inp", "")+".inp"
        inputfilepath = os.path.join(self.workfolder, inputfilename)
        myfile = open(inputfilepath, "w")
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
        if resp == "yes":
            myfile.write("resp              yes\n")
            myfile.write("esp_grid_dens     16.0\n")
            myfile.write("esp_grid_incr     0.05\n")
            myfile.write("esp_grid_layers   16  \n")

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
        return inputfilepath

    def generate_input(self, mol:System):
        charge, spinmult = 0, 1
        if isinstance(mol, Molecule):
            charge = mol.charge
        elif isinstance(mol, MoleculeComplex) or isinstance(mol, SolvatedSystem):
            charge = mol.netcharge
        spinmult = mol.multiplicity
        coord = getattr(mol, self.coord_format)
        inputfile = self.generate_tc_input(self.jobname, coord, charge, spinmult, **self.tc_args)
        self.input_file = inputfile
        
    def predict_output(self, mol:System, **kwargs):
        self.output_file = os.path.splitext(self.input_file)[0] + ".out"
        self.scr_path    = os.path.dirname(self.input_file) + "/" + "scr-" + self.jobname

    def generate_slurm_script(self, inputfile:str, **kwargs) -> str:
        # 我不知道为什么在pascal上使用srun会永远等下去但sbatch就不会。
        jobname = os.path.splitext(os.path.basename(inputfile))[0]
        outfile = jobname + ".out"
        subfile = jobname + ".sh"
        subpath = os.path.join(self.workfolder, subfile)
        myfile = open(subpath, "w")
        myfile.write("#!/bin/bash\n")
        myfile.write("#SBATCH --time=00:30:00   \n")
        myfile.write("#SBATCH --nodes=1         \n")
        myfile.write("#SBATCH --mem=1G          \n")
        myfile.write("#SBATCH --ntasks=1        \n")
        myfile.write("#SBATCH --cpus-per-task=1 \n")
        myfile.write("#SBATCH --gres=gpu:1      \n")
        myfile.write("module load {tcmodule}    \n".format(tcmodule = self.terachem_package_path))
        myfile.write("{terachem} {ifile} > {ofile}  \n".format(terachem = self.terachem, ifile = inputfile, ofile = outfile))
        myfile.close()
        return subpath
    
    def generate_cmd(self):
        if self.sbatch_use:
            self.slurm_script = self.generate_slurm_script(self.input_file, **self.tc_args)
            cmd = "sbatch {}".format(self.slurm_script)
        elif self.srun_use:
            cmd = "srun -n 1 {terachem} {ifile} > {ofile}".format(self.terachem, ifile = self.input_file, ofile = self.output_file)
        else:
            cmd = "{terachem} {ifile} > {ofile}".format(self.terachem, ifile = self.input_file, ofile = self.output_file)
            self.logger.warning("Running terachem on root node is not recommended.")
        return cmd
    
    def execute(self, cmd:str):
        if self.sbatch_use:
            j = SbatchJob(self.slurm_script, logger = self.logger)
            j.wait_until_finish()
        else:
            self.logger.info("CMD: {}".format(cmd))
            subprocess.run(cmd, shell=True, check=True)

    def check_output(self, mol:System):
        if not os.path.exists(self.output_file):
            raise RuntimeError("TeraChem output file {} not found.".format(self.output_file))

    def process_output(self, mol:System):
        pass

    # main workflow for getting result from terachem output
    def run(self, mol:System, **kwargs):
        for key, value in kwargs.items():
            self.logger.info("Set {} to {}".format(key, value))
            self.tc_args[key] = value

        self.check_system(mol)
        self.generate_input(mol, **kwargs)
        self.predict_output(mol, **kwargs)
        self.execute(self.generate_cmd())
        self.check_output(mol)
        self.process_output(mol)

    def optimize_geometry(self, mol:System, **kwargs):
        self.tc_args["run"] = "minimize"
        self.run(mol, **kwargs)
        optimxyz = os.path.join(self.scr_path, "optim.xyz")
        elems, coords = read_optimized_xyz(optimxyz)
        write_xyz(mol.xyz, elems, coords)

    def resp(self, mol:System, **kwargs):
        self.tc_args["run"] = "energy"
        self.tc_args["resp"] = "yes"
        self.run(mol, **kwargs)
        return read_resp_charge(self.output_file)
    
    def energy(self, mol:System, **kwargs):
        self.tc_args["run"] = "energy"
        self.run(mol, **kwargs)
        return read_energy(self.output_file)
    