from openbabel import openbabel as ob
from openbabel import pybel
import os, re, subprocess, shutil, glob, sys
from abc import ABC, abstractmethod
from autosolvate.globs import keywords_avail, available_qm_programs
from autosolvate.resp_classes.resp_abstract import RespABC

class RespGaussian(RespABC):
    def __init__(self, **kwargs):
        print("*"*40)
        print("Run Gaussian to generate RESP charge and mol2 file".center(40," ") )
        print("*"*40)
        self.srun_use = True
        # TODO: read in srun_use option from kwargs
        RespABC.__init__(self, **kwargs)
        
        if self.qm_exe == None:
           print("WARNING: Gaussian executable name is not specified for RESP charge fitting!")
           print("WARNING: Using g16 by default. If failed later,\
                  please rerun with the option -e specified!")
           self.qm_exe = 'g16'
        if self.qm_dir == None:
            print("WARNING: Gaussian executable directory is not specified for RESP charge fitting!")
            print("WARNING: Setting to default path: /opt/packages/gaussian/g16RevC.01/g16/")
            print("WARNING: If failed later, please rerun with the option -d specified!")
            self.qm_dir = '/opt/packages/gaussian/g16RevC.01/g16/'
            gspath = os.path.join(self.qm_dir, self.qm_exe)
            #DEBUG
            """
            if not os.path.exists(gspath):
                print("Error! Gaussian executable path",gspath)
                print("Does not exist! Exiting...")
                exit()
            """
            
    def writeGaussianInput(self):
        """
        Set up Gaussian calculation to compute electrostatic potential.
     
        optional arguments:
          calculation - one of 'optimize' (hours), or 'singlepoint' (minutes) (defaults to 'singlepoint')
     
        """
        print("-"*40)
        print("Preparing Gaussian input files...".center(40," "))
        print("-"*40)

        charge = self.molecule.GetTotalCharge()
        multiplicity = self.molecule.GetTotalSpinMultiplicity()

        gaussian_gesp = self.molname + ".gesp"
        gaussian_com = self.molname + "_gcrt.com"

        cmd1 = "$AMBERHOME/bin/antechamber -i {pdbfile} -fi pdb".format(pdbfile=self.pdbfile)
        cmd1 += " -o {com} -fo gcrt -gv 1 -ge {gesp}".format(com=gaussian_com, gesp=gaussian_gesp)
        cmd1 += " -s 2 -nc " + str(charge) + " -m " + str(multiplicity)  

        if self.srun_use:
            cmd1='srun -n 1 '+cmd1
        print(cmd1)
        subprocess.call(cmd1, shell=True)
     
    def executeGaussian(self):
        print("-"*40)
        print("Running Gaussian. This may take a while...".center(40," "))
        print("-"*40)

        gaussian_gesp = self.molname + ".gesp"
        gaussian_com = self.molname + "_gcrt.com"
        gaussian_out = self.molname + "_gcrt.out"
        
        cmd21="export GAUSS_EXEDIR=" + self.qm_dir + "; export GAUSS_SCRDIR="+self.resp_scr_dir + ";"
        #$PROJECT/TMP_GAUSSIAN;" #/expanse/lustre/projects/mit181/eh22/TMP_GAUSSIAN/;" #/scratch/$USER/$SLURM_JOBID ;"
        cmd22 = self.qm_exe + " < {com} > {out}".format(com=gaussian_com, out=gaussian_out)
        if self.srun_use:
            cmd22='srun -n 1 '+cmd22
        cmd2=cmd21+cmd22
        print(cmd2)
        subprocess.call(cmd2, shell=True)
        if not os.path.isfile(gaussian_gesp):
            print("Gaussian failed to generate solute.gesp")
            sys.stdout.flush()
            sys.exit()
        print("Gaussian ESP calculation done")
        
        

    def runGaussian(self):
        self.writeGaussianInput()
        self.executeGaussian()

    def respFit(self):
        """
        Perform RESP fit on the molecule 
        """
        print("-"*40)
        print("Start RESP fitting with gesp file generaged by Gaussian".center(40," ") )
        print("-"*40)
        
    
        # Write molecule with updated charges to mol2 file.
        gaussian_gesp = self.molname + ".gesp"
        mol2_filename = self.molname + ".mol2"

        print("Writing out the mol2 file with resp charge: " + mol2_filename)
        cmd = "$AMBERHOME/bin/antechamber"
        cmd += " -i {gesp} -fi gesp -o {mol2} -fo mol2 -c resp -eq 2 -rn SLU".format(
             gesp=gaussian_gesp, mol2=mol2_filename)

        if self.srun_use:
                 cmd='srun -n 1 '+cmd
        print(cmd)
        subprocess.call(cmd,shell=True)

    def run(self):    
        # Run Gaussian to calculate ESP potential
        if not os.path.isdir(self.resp_scr_dir):
            print("Creating the scratch folder for RESP fitting: ", self.resp_scr_dir)
            os.mkdir(self.resp_scr_dir)

        print("Copying over the pdb file", self.pdbfile, " to ", self.resp_scr_dir)
        shutil.copy(os.path.join(self.rundir,self.pdbfile), self.resp_scr_dir)

        print("Navigating to the scratch folder for RESP fitting: ", self.resp_scr_dir)
        os.chdir(self.resp_scr_dir)

        self.runGaussian()
        self.respFit()

        print("Navigating back to the folder for AutoSolvate run: ", self.rundir)
        os.chdir(self.rundir)
        
        mol2_filename = self.molname + ".mol2"
        print("Copying the generated mol2 file ", mol2_filename, " to rundir", self.rundir)
        shutil.copy(os.path.join(self.resp_scr_dir, mol2_filename), self.rundir)
