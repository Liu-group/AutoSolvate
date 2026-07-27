import autosolvate.resp_classes.gen_esp as gen_esp
from autosolvate.resp_classes.resp_abstract import RespABC
import os
import shutil
def get_crds_from_xyz(inp):
    crds = []
    with open(inp,'r') as f:
        data = f.readlines()
        for line in data[2:]:
            if len(line) > 3:
                line = line.strip()
                atom = line.split()[0].capitalize()
                x = line.split()[1]
                y = line.split()[2]
                z = line.split()[3]
                crds.append([atom,x,y,z])
    return crds


class RespORCA(RespABC):
    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.logger.info("Running ORCA to generate RESP charges and mol2 file")
        self.srun_use = kwargs.get("srun_use", False)
        
        if self.qm_exe == None:
           self.logger.warning("Orca executable not specified; defaulting to 'orca'")
           self.qm_exe = 'orca'
        if self.qm_dir == None:
            self.logger.warning("Orca directory not specified; defaulting to /opt/orca/5.0.2/")
            self.qm_dir = '/opt/orca/5.0.2/'
        
        self.orcapath = os.path.join(self.qm_dir, self.qm_exe)
        if os.path.exists(self.orcapath):
            self.logger.info(f"Found ORCA executable at {self.orcapath}")
        elif os.path.exists(self.qm_dir) and os.path.isfile(self.qm_dir) and os.access(self.qm_dir, os.X_OK):
            self.logger.info(f"Found ORCA executable at {self.qm_dir}")
            self.orcapath = self.qm_dir
        else:
            self.orcapath = self.qm_exe
            self.logger.warning(
                "ORCA executable not found at configured path; will defer check until runtime. "
                "Configured qm_dir=%s qm_exe=%s",
                self.qm_dir,
                self.qm_exe,
            )

    def _require_executable(self, exe_name: str):
        if os.path.isabs(exe_name):
            if os.path.exists(exe_name) and os.access(exe_name, os.X_OK):
                return exe_name
        else:
            resolved = shutil.which(exe_name)
            if resolved is not None:
                return resolved
        raise FileNotFoundError(f"Required executable not found at runtime: {exe_name}")

    def writeOrcaInput(self):
        """
        Set up Orca calculation to compute electrostatic potential.
     
        optional arguments:
          calculation - one of 'optimize' (hours), or 'singlepoint' (minutes) (defaults to 'singlepoint')
     
        """
        self.logger.info("Preparing Orca input files")

        charge = self.molecule.GetTotalCharge()
        multiplicity = self.molecule.GetTotalSpinMultiplicity()

        orca_inp = self.molname + "_orca.in"
        orca_out = self.molname + '_orca.out'
        xyzfile = self.xyzfile
        orca_crd = get_crds_from_xyz(xyzfile)
        
        with open(orca_inp,'w') as ofile:
            ofile.write('! hf 6-31g* keepdens\n\n')
            ofile.write("%pal\nnprocs " + str(self.nprocs) + "\nend\n\n")
            ofile.write('* xyz ' + str(charge) + ' ' + str(multiplicity) + '\n')
            for line in orca_crd:
                atomname = line[0]
                x  = line[1]
                y = line[2]
                z = line[3]
                ofile.write(atomname + '  ' + x + '  ' + y  + '  ' + z + '\n')
            ofile.write('*')
            
        cmd1 = self.orcapath + ' ' + orca_inp + ' > ' + orca_out
        
        ## check convergence of orcaout ##
        convergence = False
        if os.path.exists(orca_out):
            with open(orca_out) as f:
                data = f.read()
                if 'ORCA TERMINATED NORMALLY' in data:
                    self.logger.info('Orca calculation finished for ESP generation')
                    convergence = True
            
        if not convergence:
            self._require_executable(self.orcapath)
            if self.srun_use:
                procs = str(self.nprocs or 1)
                cmd1='srun -n ' + procs + ' ' +cmd1
            self.run_shell_command(
                cmd1,
                "Running ORCA ESP calculation",
                produced_file=os.path.abspath(orca_out),
                purpose="ORCA output log",
            )
    
    def generateESP(self):
        vpotpath = self.orcapath + '_vpot'
        orca_inp = self.molname + "_orca.in"
        orca_out = self.molname + '_orca.out'
        gbw = self.molname + '_orca.gbw'
        density = self.molname + '_orca.scfp'
        espf =  self.molname + '_orca.esp'
        out =  self.molname + '_orca.espout'

        if os.path.exists(espf) and os.path.getsize(espf) > 0:
            self.logger.info("Found existing ESP file %s; skipping ORCA vpot generation", os.path.abspath(espf))
            return

        self._require_executable(vpotpath)

        crds = get_crds_from_xyz(self.xyzfile)
        elements = []
        newcrds = []
        for line in crds:
            newcrds.append([float(line[1]),float(line[2]),float(line[3])])
            elements.append(line[0].capitalize())
        
        gen_esp.gen_grids(crds=newcrds,elements=elements,orcapath=vpotpath,gbw=gbw,denisty=density,out=out)
        gen_esp.convertoesp(espin=out,espout=espf,crds=newcrds)
        self.logger.info("The ORCA ESP grid is generated at %s", os.path.abspath(espf))
        
       # os.system("resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg \
       #        -e %s -s resp1_calc.esp" %espf)
      #  os.system("resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg \
       #       -t resp2.chg -e %s -s resp2_calc.esp" %espf)
        
    
    
    def respFit(self):
        Method1 = True
        self.logger.info("Start RESP fitting with Orca ESP output")
        
        ac_filename = self.molname + ".ac"
        self.logger.info("Generating Antechamber file: %s", ac_filename)
        
        cmd = (
            "antechamber -i {pdbfile} -fi pdb -o {ac} -fo ac -c dc -nc {charge} -pl 30 -rn {resname}".format(
                pdbfile=self.pdbfile,
                ac=ac_filename,
                charge=self.charge,
                resname=self.residue_name,
            )
        )
        self.run_shell_command(
            cmd,
            "Generating Antechamber file",
            log_file=os.path.join(self.rundir, 'antechamber_acfile.log'),
            produced_file=os.path.abspath(ac_filename),
            purpose="Antechamber intermediate",
        )
                           
        resp1_filename = self.molname + '.resp1'
        self.logger.info("Generating stage 1 RESP input: %s", resp1_filename)
        
        cmd = "respgen -i {ac} -o {resp1} -f resp1 -l 10".format(ac=ac_filename, resp1=resp1_filename)
        self.run_shell_command(
            cmd,
            "Generating stage 1 RESP input",
            produced_file=os.path.abspath(resp1_filename),
            purpose="Stage 1 RESP input deck",
        )

        # Stage 2
        resp2_filename = self.molname + '.resp2'
        self.logger.info("Generating stage 2 RESP input: %s", resp2_filename)
        
        cmd = "respgen -i {ac} -o {resp2} -f resp2 -l 10".format(ac=ac_filename, resp2=resp2_filename)
        self.run_shell_command(
            cmd,
            "Generating stage 2 RESP input",
            produced_file=os.path.abspath(resp2_filename),
            purpose="Stage 2 RESP input deck",
        )
        
        esp_filename =  self.molname + '_orca.esp'
        
        # Perform RESP fit.
        #Stage 1
        qout1_filename = self.molname + '.qout_stage1'
        respout1_filename = self.molname + '.respout1'
        self.logger.info("Running stage 1 RESP fit -> %s", qout1_filename)
        
        cmd = "resp -O -i {resp1} -o {respout1} -e {esp} -t {qout1}".format(
            resp1=resp1_filename, respout1=respout1_filename,
            esp=esp_filename, qout1=qout1_filename)
        self.run_shell_command(
            cmd,
            "Running stage 1 RESP fitting",
            produced_file=os.path.abspath(qout1_filename),
            purpose="Stage 1 RESP charges",
        )

        #Stage 2
        qout2_filename = self.molname + '.qout_stage2'
        respout2_filename = self.molname + '.respout2'
        self.logger.info("Running stage 2 RESP fit -> %s", qout2_filename)

        cmd = "resp -O -i {resp2} -o {respout2} -e {esp} -q {qout1} -t {qout2}".format(
            resp2=resp2_filename, respout2=respout2_filename, esp=esp_filename,
            qout1=qout1_filename, qout2=qout2_filename)
        self.run_shell_command(
            cmd,
            "Running stage 2 RESP fitting",
            produced_file=os.path.abspath(qout2_filename),
            purpose="Stage 2 RESP charges",
        )

        # Write molecule with updated charges to mol2 file.
        mol2_filename = self.molname + ".mol2"
        self.logger.info("Writing mol2 with RESP charges: %s", mol2_filename)
        cmd = (
            "antechamber -i {ac} -fi ac -o {mol2} -fo mol2 -c rc -cf {qout2} -pl 30 -rn {resname}".format(
                ac=ac_filename,
                mol2=mol2_filename,
                qout2=qout2_filename,
                resname=self.residue_name,
            )
        )
        self.run_shell_command(
            cmd,
            "Finalizing RESP mol2 output",
            produced_file=os.path.abspath(mol2_filename),
            purpose="RESP mol2 file",
        )
        
 
         #   with open('resp1.in','w') as fresp1:
          #      fresp1.write()
            
    def run(self):
      #  if not os.path.isdir(self.resp_scr_dir):
       #     print("Creating the scratch folder for RESP fitting: ", self.resp_scr_dir)
       #     os.mkdir(self.resp_scr_dir)

      #  print("Copying over the pdb file", self.pdbfile, " to ", self.resp_scr_dir)
      #  shutil.copy(os.path.join(self.rundir,self.pdbfile), self.resp_scr_dir)

      #  print("Navigating to the scratch folder for RESP fitting: ", self.resp_scr_dir)
      #  os.chdir(self.resp_scr_dir)
        
        self.writeOrcaInput()
        self.generateESP()
        self.respFit()

     #   print("Navigating back to the folder for AutoSolvate run: ", self.rundir)
     #   os.chdir(self.rundir)
        
     #   mol2_filename = self.molname + ".mol2"
    #    print("Copying the generated mol2 file ", mol2_filename, " to rundir", self.rundir)
     #   shutil.copy(os.path.join(self.resp_scr_dir, mol2_filename), self.rundir)
