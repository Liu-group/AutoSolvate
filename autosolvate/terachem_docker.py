#--------------------------------------------------------------------------------------------------#
# terachem_docker.py. 
# Description: 
#   This module calls TeraChem for RESP charge fitting. 
#   This is just to avoid the situation where Gaussian is not installed.
#--------------------------------------------------------------------------------------------------#
import os, subprocess
import time
from typing import List

class SbatchJob():
    def __init__(self, script:str, maxtime = 86400*7, submit_dir = "", check_interaval = 10):
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
            print("warning: the time " + s + " cannot be processed.")
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
            print(f"failed to submit script {self.name}. slurm report:")
            print(result)

    def wait_until_finish(self):
        if self.status == "U":
            self.submit()
        for i in range(self.maximum_time // self.check_interval):
            time.sleep(self.check_interval)
            self.check_finished()
            if self.status == "F":
                print(f"job {self.jobid:9d} {self.jobname:9s} finished. ")
                break
            else:
                print(f"job {self.jobid:9d} {self.jobname:9s} status {self.status:2s} time {self.time:9d}")
        else:
            print(f"job {self.jobid:9d} {self.jobname:9s} has been running for more than {self.maximum_time} seconds, no longer waiting for this task to complete.")


class TeraChemDocker(object):
    """
    A temporary function class that supports RESP charge fitting using terachem. In case Gaussuian is not install but TeraChem is.
    """

    def __init__(self, terachem_path:str="") -> None: 
        '''
        @NOTE: 
        Formatted as Patrick's code.

        '''
        self.terachem = terachem_path
        self.terachem_package_path = ""
        self.get_tc_path(terachem_path)

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
                        raise ModuleNotFoundError("No terachem in: " + " ".join(traversed_path))
        self.terachem = terachem_path

    def generate_tc_input(self, xyzname:str, charge:int, spinmult:int, jobetype:str = "resp", precision:str = "dynamic") -> str:
        """
        pass
        """
        name = os.path.splitext(os.path.basename(xyzname))[0]

        jobname = name + "-respfit"
        myfile = open(jobname+".inp", "w")
        # necessary information
        myfile.write("jobname        " + jobname + "\n")
        myfile.write("coordinates    " + xyzname + "\n")
        myfile.write("run            energy\n")
        myfile.write("scrdir         scr-" + jobname + "\n")
        myfile.write("gpus           1\n")

        # necessary information for calculation
        myfile.write("method         uhf\n")
        myfile.write("basis          6-31gs\n") 
        myfile.write("charge         " + str(charge) + "\n")
        myfile.write("spinmult       " + str(spinmult) + "\n")
        myfile.write("sphericalbasis yes\n")

        # scf setting
        myfile.write("guess          generate\n")
        myfile.write("maxit          1000\n")
        myfile.write("scf            diis+a\n")
        myfile.write("precision      " + str(precision) + "\n")
        myfile.write("convthre       3.0e-5\n")
        myfile.write("xtol           1.0e-6\n")

        # specified setting
        if jobetype == "resp":
            myfile.write("resp              yes\n")
            myfile.write("esp_grid_dens     16.0\n")
            myfile.write("esp_grid_incr     0.05\n")
            myfile.write("esp_grid_layers   16  \n")
        elif jobetype == "minimize":
            myfile.write("min_coordinates   tric\n")
            myfile.write("new_minimizer     yes\n")
            myfile.write("min_converge_e    1.0e-6\n")
            myfile.write("purify            no\n")
            myfile.write("nstep             1000\n")

        myfile.write("\nend\n")
        myfile.close()

        return jobname+".inp"

    def generate_srun_cmd(self, inputfile:str) -> str:
        cmd = "srun -n 1 {terachem} {ifile}.inp > {ifile}.out".format(self.terachem, ifile = inputfile)
        return cmd
    
    def run_cmd(self, cmd:str) -> None:
        subprocess.call(cmd, shell = True)

    def generate_slurm_script(self, inputfile:str) -> str:
        # 我不知道为什么在pascal上使用srun会永远等下去但sbatch就不会。
        jobname = os.path.splitext(os.path.basename(inputfile))[0]
        outfile = jobname + ".out"
        subfile = jobname + ".sh"
        myfile = open(subfile, "w")
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
        return subfile
    
    def run_sbatch(self, slurmscript:str) -> None:
        # only used when submit a slurm script instead of using srun
        sjob = SbatchJob(script=slurmscript, maxtime=60*60, check_interaval=5)
        sjob.submit()
        sjob.wait_until_finish()

    def get_resp_charge(self, outfile:str) -> list:
        with open(outfile, "r") as f:
            content = f.read()
        lines = content.splitlines()
        for i, line in enumerate(lines):
            if line.find("ESP restraint charges:") != -1:
                break
        newlines = []
        for line in lines[i+3:]:
            if line.startswith("-----------------"):
                break
            newlines.append(line)
        chgs = []
        for line in newlines:
            data = line.split()
            chg = float(data[4])
            chgs.append(chg)
        return chgs

    def get_optim_geometry(self, trjxyz:str) -> None:
        '''
        @ISSUE:
        If the input file is in pdb format, TeraChem will also write a pdb file with its optimized structure. However, there are data errors in the pdb file, so it cannot be read correctly.
        '''
        pass

class Mol2Faker(object):
    """
    Antechamber is not compatible with terachem. This class will first adjust the charge of the molecule and make it possible to use the bcc method for charge calculation, then call Antechamber to generate a template mol2 file, and then manually modify the atomic charges in the file.
    """
    def __init__(self, pdbname:str, resname:str, charge:int, spinmult:str):
        self.pdbname    = pdbname
        self.name       = os.path.splitext(pdbname)[0]
        self.resname    = resname
        self.mol2name   = self.name + ".mol2"
        self.fakemol2   = self.name + "-fake.mol2"
        self.charge     = charge
        self.spinmult   = spinmult
        self.fakecharge = charge
        self.fakemult   = 1

    def removeConectFromPDB(self) -> None:
        print("cleaning up {pdbname}".format(pdbname = self.pdbname))
        pdb1 = open(self.pdbname).readlines()
        pdb2 = open(self.pdbname,'w')
        for line in pdb1:
            if 'CONECT' not in line:
                pdb2.write(line)
        pdb2.close()

    def calculate_fake_charge(self) -> int:
        """create a fake charge that can make the spinmult to be 1."""
        if self.spinmult <= 1:
            return self.charge
        if self.spinmult % 2 == 1:
            return self.charge
        if self.charge <= 0:
            self.fakecharge = self.charge + 1
        elif self.charge > 0:
            self.fakecharge = self.charge - 1

    def fake_mol2(self) -> None:
        self.removeConectFromPDB()
        self.calculate_fake_charge()
        pl = -1
        cmd5 = f"$AMBERHOME/bin/antechamber -i {self.pdbname} -fi pdb -o {self.fakemol2} -fo mol2 -c bcc -eq 2 -rn {self.resname} -pl {pl} -nc {self.fakecharge} -m {self.fakemult}"
        subprocess.call(cmd5, shell = True)

    def write_mol2_line(self, res) -> str:
        newline = ""
        newline += " {:>3d}"     .format(int(res[0]))
        newline += " {:<4s}"     .format(res[1])
        newline += " {:>12.6f} {:>12.6f} {:>12.6f}".format(float(res[2]), float(res[3]), float(res[4]))
        newline += " {:<6s}"    .format(res[5])
        newline += " {:>3d}"     .format(int(res[6]))
        newline += " {:<7s}"    .format(res[7])
        newline += " {:>9.6f} "  .format(float(res[8]))
        if len(res) > 9:
            newline += " ".join(res[9:])
        newline += "\n"
        return newline

    def modify_mol2(self, charges:List[float]) -> None:
        f1 = open(self.fakemol2, "r")
        f2 = open(self.mol2name, "w")
        section = ""
        atomid = 0
        for line in f1:
            if line.startswith("@"):
                section = line.strip("\n").replace("@<TRIPOS>", "")
            if section.find("ATOM") == -1:
                f2.write(line)
            else:
                res = line.split()
                if len(res) == 9 or len(res) == 10:
                    res[8] = str(charges[atomid])
                    newline = self.write_mol2_line(res)
                    f2.write(newline)
                    atomid += 1
                else:
                    f2.write(line)
        f1.close()
        f2.close()

def check_terachem(tcpath:str = "terachem"):
    tcdoker     = TeraChemDocker(tcpath)

def resp_terachem(pdbname:str, resname:str, charge:int, spinmult:int, tcpath:str = "") -> None:
    tcdoker     = TeraChemDocker(tcpath)
    inputfile   = tcdoker.generate_tc_input(pdbname, charge, spinmult)
    slurmscript = tcdoker.generate_slurm_script(inputfile)
    tcdoker.run_sbatch(slurmscript)
    chgs        = tcdoker.get_resp_charge(os.path.splitext(inputfile)[0] + ".out")
    
    m2faker     = Mol2Faker(pdbname, resname, charge, spinmult)
    m2faker.fake_mol2()
    m2faker.modify_mol2(charges = chgs)

if __name__ == "__main__":
    resp_terachem("benzene.pdb", "BEZ", 1, 2)

