import os
import shutil
import sys

from autosolvate.resp_classes.resp_abstract import RespABC

class RespGaussian(RespABC):
    def __init__(self, **kwargs):
        RespABC.__init__(self, **kwargs)
        self.logger.info("Running Gaussian to generate RESP charge and mol2 file")
        self.srun_use = kwargs.get("srun_use", True)
        
        if self.qm_exe is None:
            self.logger.warning("Gaussian executable not specified; defaulting to g16")
            self.qm_exe = 'g16'
        if self.qm_dir is None:
            self.qm_dir = '/opt/packages/gaussian/g16RevC.01/g16/'
            self.logger.warning("Gaussian directory not specified; defaulting to %s", self.qm_dir)
            
    def writeGaussianInput(self):
        """
        Set up Gaussian calculation to compute electrostatic potential.
     
        optional arguments:
          calculation - one of 'optimize' (hours), or 'singlepoint' (minutes) (defaults to 'singlepoint')
     
        """
        self.logger.info("Preparing Gaussian input files")

        charge = self.molecule.GetTotalCharge()
        multiplicity = self.molecule.GetTotalSpinMultiplicity()

        gaussian_gesp = self.molname + ".gesp"
        gaussian_com = self.molname + "_gcrt.com"

        cmd1 = "$AMBERHOME/bin/antechamber -i {pdbfile} -fi pdb".format(pdbfile=self.pdbfile)
        cmd1 += " -o {com} -fo gcrt -gv 1 -ge {gesp}".format(com=gaussian_com, gesp=gaussian_gesp)
        cmd1 += " -s 2 -nc " + str(charge) + " -m " + str(multiplicity)
        cmd1 += f" -rn {self.residue_name}"

        if self.srun_use:
            cmd1='srun -n 1 '+cmd1
        self.run_shell_command(
            cmd1,
            "Generating Gaussian input via Antechamber",
            produced_file=os.path.abspath(gaussian_com),
            purpose="Gaussian input deck",
        )
     
    def executeGaussian(self):
        self.logger.info("Running Gaussian job")

        gaussian_gesp = self.molname + ".gesp"
        gaussian_com = self.molname + "_gcrt.com"
        gaussian_out = self.molname + "_gcrt.out"
        
        cmd21="export GAUSS_EXEDIR=" + self.qm_dir + "; export GAUSS_SCRDIR="+self.resp_scr_dir + ";"
        #$PROJECT/TMP_GAUSSIAN;" #/expanse/lustre/projects/mit181/eh22/TMP_GAUSSIAN/;" #/scratch/$USER/$SLURM_JOBID ;"
        cmd22 = self.qm_exe + " < {com} > {out}".format(com=gaussian_com, out=gaussian_out)
        if self.srun_use:
            cmd22='srun -n 1 '+cmd22
        cmd2=cmd21+cmd22
        self.run_shell_command(
            cmd2,
            "Executing Gaussian ESP job",
            produced_file=os.path.abspath(gaussian_out),
            purpose="Gaussian output log",
        )
        if not os.path.isfile(gaussian_gesp):
            raise RuntimeError("Gaussian failed to generate .gesp file")
        self.logger.info("The Gaussian ESP file is generated at %s", os.path.abspath(gaussian_gesp))
        
        

    def runGaussian(self):
        self.writeGaussianInput()
        self.executeGaussian()

    def respFit(self):
        """
        Perform RESP fit on the molecule 
        """
        self.logger.info("Start RESP fitting with Gaussian ESP output")
        
    
        # Write molecule with updated charges to mol2 file.
        gaussian_gesp = self.molname + ".gesp"
        mol2_filename = self.molname + ".mol2"

        self.logger.info("Writing mol2 file with RESP charges: %s", mol2_filename)
        cmd = "$AMBERHOME/bin/antechamber"
        cmd += " -i {gesp} -fi gesp -o {mol2} -fo mol2 -c resp -eq 2 -rn {resname}".format(
            gesp=gaussian_gesp, mol2=mol2_filename, resname=self.residue_name)

        if self.srun_use:
            cmd='srun -n 1 '+cmd
        self.run_shell_command(
            cmd,
            "Converting Gaussian ESP to RESP mol2",
            produced_file=os.path.abspath(os.path.join(self.resp_scr_dir, mol2_filename)),
            purpose="RESP mol2 file",
        )

    def run(self):    
        # Run Gaussian to calculate ESP potential
        os.makedirs(self.resp_scr_dir, exist_ok=True)
        self.logger.info("Using RESP scratch directory: %s", self.resp_scr_dir)
        shutil.copy(os.path.join(self.rundir,self.pdbfile), self.resp_scr_dir)
        self.logger.debug("Copied %s to scratch", self.pdbfile)

        cwd = os.getcwd()
        try:
            os.chdir(self.resp_scr_dir)
            self.runGaussian()
            self.respFit()
        finally:
            os.chdir(cwd)
        
        mol2_filename = self.molname + ".mol2"
        shutil.copy(os.path.join(self.resp_scr_dir, mol2_filename), self.rundir)
        self.logger.info("RESP mol2 copied to %s", self.rundir)