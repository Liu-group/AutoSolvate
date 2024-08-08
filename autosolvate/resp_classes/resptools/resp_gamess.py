from openbabel import openbabel as ob
from openbabel import pybel
import os, re, subprocess, shutil, glob, sys
from abc import ABC, abstractmethod
from autosolvate.resp_classes.resp_abstract import RespABC


class RespGAMESS(RespABC):
    def __init__(self, **kwargs):
        print("*"*40)
        print("Run GAMESS to generate RESP charge and mol2 file".center(40," ") )
        print("*"*40)
        self.srun_use = kwargs["srun_use"] if "srun_use" in kwargs else True # GAMESS script usually requires srun_use
        self.nnodes = kwargs["nnodes"] if "nnodes" in kwargs else 1 # Gamess need
        self.ncpus = kwargs["ncpus"] if "ncpus" in kwargs else 1
        self.version = kwargs["gamessversion"] if "gamessversion" in kwargs else "01"
        RespABC.__init__(self, **kwargs)
        self.potential = None

        if self.qm_exe == None:
           print("WARNING: Gamess executable name is not specified for RESP charge fitting!")
           print("WARNING: Using rungms by default. If failed later,\
                  please rerun with the option -e specified!")
           self.qm_exe = 'rungms'
        
        if self.qm_dir == None:
           print("WARNING: GAMESS executable directory is not specified for RESP charge fitting!")
           print("WARNING: Setting to default path: ")
           print("WARNING: If failed later, please rerun with the option -d specified!")
           self.qm_dir = '/opt/gamess/'
           gamess_path = os.path.join(self.qm_dir, self.qm_exe)
           if not os.path.exists(gamess_path):
               print("Error! Gamess executable path",gamess_path)
               print("Does not exist! Exiting...")
               exit()
            

    def findLineStartWith(self, lines, string):
        """\
        Find the first line that begins with the given string.
    
        required arguments:
          lines - the lines to search
          string - the given string to check for
    
        returns:
          index - the index of the line (starting from zero) where the match was found, 
	          or 'None' if no match found
        """
    
        index = 0
        for line in lines:
            if line.startswith(string):
                return index
            index += 1
    
        raise Exception("Line not found.")   

    def findLineContains(self, lines, string):
        """\
        Find the first line that contains the given string.
    
        required arguments:
          lines - the lines to search
          string - the given string to check for
    
        returns:
          index - the index of the line (starting from zero) where the match was found, 
	          or 'None' if no match found
        """
    
        index = 0
        for line in lines:
            if string in line:
                return index
            index += 1
    
        raise Exception("Line not found.")   

    def readFile(self, filename):
        """
        Read the text contents of the specified file.
    
        required arguments:
          filename - the filename of the file whose contents are to be read
    
        return values:
          lines - the contents of the file, as a list of lines.
        """
        
        # Open the file for reading in text mode.
        infile = open(filename, 'rt')
    
        # Read contents of file to a list.
        lines = infile.readlines()
          
        # Return lines.
        return lines

    def writeGAMESSInput(self, calculation = 'singlepoint'):
        """
        Set up GAMESS calculation to compute electrostatic potential.
     
        optional arguments:
          calculation - one of 'optimize' (hours), or 'singlepoint' (minutes) (defaults to 'singlepoint')
     
        """
        print("-"*40)
        print("Preparing GAMESS input files...".center(40," "))
        print("-"*40)

        charge = self.molecule.GetTotalCharge()
        multiplicity = self.molecule.GetTotalSpinMultiplicity()
        scf = "RHF" if multiplicity==1 else "UHF"
     
        # GAMESS header blocks.
        # GAMESS is instructed to carry out a RHF or UHF claculation with a 6-31G* basis set, after which
        # the electrostatic potential on concentric Connolly surfaces is evaluated.
        gamess_header = { }
        
        # Header for energy minimization.
        gamess_header['optimize'] = """\
 $CONTRL SCFTYP={scf} EXETYP=RUN RUNTYP=OPTIMIZE COORD=UNIQUE $END
 $CONTRL ICHARG={charge} MULT={multiplicity} $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS NPRINT=8 $END
 $SCF    DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $STATPT NSTEP=100 NPRT=-2 NPUN=-2 $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
""".format(scf=scf, charge=charge, multiplicity=multiplicity)
     
        # Header for no optimization case.
        gamess_header['singlepoint'] = """\
 $CONTRL SCFTYP={scf} EXETYP=RUN RUNTYP=ENERGY COORD=UNIQUE $END
 $CONTRL ICHARG={charge} MULT={multiplicity} $END
 $CONTRL MOLPLT=.TRUE. UNITS=ANGS NPRINT=8 $END
 $SCF    DIRSCF=.TRUE. FDIFF=.FALSE. $END
 $BASIS  GBASIS=N31 NGAUSS=6 NDFUNC=1 $END
 $ELPOT  IEPOT=1 WHERE=PDC OUTPUT=PUNCH $END
 $PDC    PTSEL=CONNOLLY CONSTR=NONE $END
 $GUESS  GUESS=HUCKEL $END
""".format(scf=scf, charge=charge, multiplicity=multiplicity)
     
     
        # Check calculation option.
        if calculation not in gamess_header.keys():
            raise ParameterError("Optional argument 'calculation' must be one of: " + gamess_header.keys())
        
        # Construct the GAMESS input file.
        # Select header for appropriate calculation.
        gamess_input = gamess_header[calculation]
        # Append data block containing Cartesian atomic coordinates.
        gamess_input += " $DATA\n"  # block header
        gamess_input += "solute\n"  # molecule name        
        gamess_input += "C1\n" # symmetry point group        
        for i in range(1,self.molecule.NumAtoms()+1):  # atomic coordinates
            atom = self.molecule.GetAtom(i)
            #elem = re.sub(r'[0-9]+', '', atom.GetType())
            elem = atom.GetType()
            gamess_input += "%10s %4.0f %13f %13f %13f\n" % (elem, atom.GetAtomicNum(), 
                             atom.GetX(), atom.GetY(), atom.GetZ())
        gamess_input += " $END\n"  # terminate block
     
        # Write input file to work directory.
        gamess_input_path = self.molname + "_gamess.inp"
        finput = open(gamess_input_path,'w')
        finput.write(gamess_input)
        finput.close()

    def readGAMESSOutput(self):
        """\
        Process GAMESS output to extract optimized coordinates, Mulliken charges, and electrostatic potential.
    
        """
    
        print("Reading GAMESS output from %s..." % self.resp_scr_dir)
        
        # Read in file contents.
        log_lines = self.readFile(os.path.join(os.path.abspath(self.resp_scr_dir), self.molname + '_gamess.log'))
        dat_lines = self.readFile(os.path.join(os.path.abspath(self.resp_scr_dir),self.molname + '_gamess.dat'))
    
        # Determine number of atoms.
        natoms = self.molecule.NumAtoms()
    
        # Extract electrostatic potential.
        import re
        index = self.findLineStartWith(log_lines, ' NUMBER OF POINTS SELECTED FOR FITTING =')
        npoints = int(re.match(' NUMBER OF POINTS SELECTED FOR FITTING =\s*(\d+)\n', log_lines[index]).groups()[0])
        index = self.findLineStartWith(dat_lines, ' ELECTROSTATIC POTENTIAL, IPT,X,Y,Z,ELPOTT')
        self.potential = { } # Electrostatic potential: potential[(x,y,z)] = psi
        for line in dat_lines[index+2:index+npoints+2]:            
            # Parse line to extract coordinate triple (x,y,z) and potential psi in bohr units.
            # Fortran format: (1X,I5,2X,3F10.5,E15.6)
            # Example format: 
            #     1     1.59920   3.28673   1.32987  -0.990640E-02
            x = float(line[8:18])
            y = float(line[18:28])
            z = float(line[28:38])
            psi = float(line[38:54])
            
            # Store potential at this grid point.
            self.potential[(x,y,z)] = psi
    
        # Check number of points read.
        if (len(self.potential.keys()) != npoints):
            raise Exception( "Only read %d points (expected %d)." % (len(self.potential.keys()), npoints) )
            
        # Extract Mulliken charges.
        # Example format:
        #          1         2         3         4         5
        #012345678901234567890123456789012345678901234567890
        #C1           6.16547  -0.16547   6.10684  -0.10684        
        index = self.findLineStartWith(dat_lines, ' POPULATION ANALYSIS')
        mulliken_charges = {}  # Mulliken charges
        idx=1
        for line in dat_lines[index+1:index+1+natoms]:
            # Parse line.
            mulliken_population = float(line[10:20])
            mulliken_charge = float(line[20:30])
            lowdin_population = float(line[30:40])
            lowdin_charge = float(line[40:50])
                                                  
            # Store Mulliken charge.
            mulliken_charges[idx] = mulliken_charge
            idx += 1
        
        # Assign charges to atoms in the OBMol object
        for i in range(1, self.molecule.NumAtoms()+1):
            atom = self.molecule.GetAtom(i)
            atom.SetPartialCharge(mulliken_charges[i])
    
        # Extract energy.
        # Example format:
    #0123456789012345678901234567890123456789012345678901234567890123456789
    # FINAL RHF ENERGY IS      -78.0308933274 AFTER   8 ITERATIONS
        pattern = ' FINAL RHF ENERGY IS'
        if self.molecule.GetTotalSpinMultiplicity() > 1:
            pattern =' FINAL UHF ENERGY IS'
        index = self.findLineStartWith(log_lines, pattern)
        hartree_to_kcal_per_mole = 627.509391
        line = log_lines[index]
        energy = float(line[20:39]) * hartree_to_kcal_per_mole
        self.molecule.SetEnergy(energy)
        

    def writeESPFile(self, esp_filename):
        
        # Write potential in RESP format.
        espfile = open(esp_filename, 'w')
        
        print("Reading optimized geometries and electrostatic potentials from GAMESS output and writing potential file...")
        # Read GAMESS output files to obtain final optimized geometry and electrostatic potential.
        self.readGAMESSOutput()
        # Write to ESP file.
        npoints = len(self.potential.keys())
        espfile.write('%5d%5d\n' % (self.molecule.NumAtoms(), npoints))
        # Write atom coordinates in Fortran (17x,3e16.7) format (in bohrs).
        for i in range(1,self.molecule.NumAtoms()+1):
            atom = self.molecule.GetAtom(i)
            # Convert from Angstroms to Bohr
            angstroms_per_bohr = 0.52917725            
            point = ()
            point += ( atom.GetX() / angstroms_per_bohr, )
            point += ( atom.GetY() / angstroms_per_bohr, )
            point += ( atom.GetZ() / angstroms_per_bohr, )
            # Write
            espfile.write('%17s%16.7e%16.7e%16.7e\n' % (('',) + point) )
        # Write potential in Fortran (qpot,x,y,z) in (1x,4e16.7) format
        for point in self.potential.keys():
            espfile.write(' %16.7e%16.7e%16.7e%16.7e\n' % ((self.potential[point],) + point) )
        espfile.close()

    def executeGAMESS(self):
        print("-"*40)
        print("Running GAMESS. This may take a while...".center(40," "))
        print("-"*40)
        if not "GAMESS_VERSION" in os.environ:
            gamess_version=2022
        else:
            gamess_version = int(os.environ["GAMESS_VERSION"])

        # Note down the current directory
        current_dir = os.getcwd()
        gamess_log = os.path.abspath(os.path.join(os.path.abspath(self.resp_scr_dir), self.molname + "_gamess.log"))

        if gamess_version<2022:
            self.version='00'
            gamess_inp = self.molname+"_gamess.inp"
            gamess_dat = self.molname + "_gamess.dat"
            os.chdir(os.path.abspath(self.resp_scr_dir))
            cmd = os.path.join(self.qm_dir, self.qm_exe) + " " \
                  + gamess_inp + " " \
                  + self.version + " " \
                  + str(self.nnodes * self.ncpus) + " " \
                  + " > " + gamess_log
        else:
            gamess_inp = os.path.abspath(os.path.join(os.path.abspath(self.resp_scr_dir),self.molname + "_gamess.inp"))



            cmd = os.path.join(self.qm_dir, self.qm_exe) + " " \
                + gamess_inp + " " \
                + self.version + " " \
                + str(self.nnodes * self.ncpus) + " " \
                + str(self.ncpus) + " " \
                + "0" + " " \
                + "0" + " " \
                + str(os.path.abspath(self.resp_scr_dir)) + " " \
                + str(os.path.abspath(self.resp_scr_dir)) + " " \
                + str(os.path.abspath(self.qm_dir)) \
                + " > " + gamess_log

        if self.srun_use:
            cmd='srun -n 1 '+cmd

        print(cmd)
        subprocess.call(cmd, shell=True)
        # Bring back the current working directory to earlier directory before the GAMESS cmd was run
        os.chdir(current_dir)
        if not os.path.isfile(gamess_log):
            print("GAMESS failed to generate log file")
            sys.stdout.flush()
            sys.exit()


        print("Copy back GAMESS .dat file from SCR")
        gamess_log_content = self.readFile(gamess_log)
        
        try:
             lineid = self.findLineStartWith(gamess_log_content, " ddikick.x: exited gracefully.")
        except:
             raise Exception("GAMESS ESP calcualtion failed."
                             + f" Please check the log file: {gamess_log}"
                             + " for details")
        if gamess_version<2022:
            dat_path = os.path.join(os.path.abspath(self.resp_scr_dir),"scr",gamess_dat)
            
        else:
            line_dat = lineid + self.findLineContains(gamess_log_content[lineid:], gamess_dat)
            dat_path = gamess_log_content[line_dat].split(" ")[-1].strip('\n')
        shutil.copy(dat_path,os.path.abspath(self.resp_scr_dir))
        line_node = self.findLineStartWith(gamess_log_content, "This job is running on host")
        node_name = gamess_log_content[line_node].split(" ")[-1].strip('\n')
        if self.srun_use:
            cmd = "scp " + node_name + ":" + dat_path + " " + os.path.abspath(self.resp_scr_dir)
            print(cmd)
            subprocess.call(cmd, shell=True)
        
        

    def runGAMESS(self):
        self.writeGAMESSInput()
        self.executeGAMESS()

    def respFit(self):
        """
        Perform RESP fit on the molecule 
    
        required arguments:
          directory - the directory to process
    
        optional arguments:
          net_charge - integral net charge of the molecule
    
        """
        print("-"*40)
        print("Start RESP fitting with ESP generaged by GAMESS".center(40," ") )
        print("-"*40)
        
    
        # Create an Antechamber file from the pdb file without assigning charges (all zeros)
        ac_filename = self.molname + ".ac"
        print("Generating Antechamber file: " + ac_filename)

        cmd = "antechamber -i {pdbfile} -fi pdb -o {ac} -fo ac -c dc -nc {charge}".format(pdbfile=self.pdbfile, ac=ac_filename, charge=self.charge)
        print(cmd)
        subprocess.call(cmd,shell=True)
        # TODO handle exceptions
        
        # Use Antechamber's 'respgen' to generate RESP input files for stage 1 and stage 2 fitting.
        # Stage 1
        resp1_filename = self.molname + '.resp1'
        print("Using respgen to generate stage 1 RESP fitting input file: " + resp1_filename)

        cmd = "respgen -i {ac} -o {resp1} -f resp1".format(ac=ac_filename, resp1=resp1_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)
        resp1_contents = self.readFile(resp1_filename)
        
        # Stage 2
        resp2_filename = self.molname + '.resp2'
        print("Using respgen to generate stage 2 RESP fitting input file: " + resp2_filename)

        cmd = "respgen -i {ac} -o {resp2} -f resp2".format(ac=ac_filename, resp2=resp2_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)
        resp2_contents = self.readFile(resp2_filename)
    
    
        # Write potential in RESP format.
        esp_filename = self.molname + ".esp"
        self.writeESPFile(esp_filename)
    
        # Perform RESP fit.
        #Stage 1
        qout1_filename = self.molname + '.qout_stage1'
        respout1_filename = self.molname + '.respout1'
        print("Using resp to perform stage 1 RESP fitting to generate: " + qout1_filename)

        cmd = "resp -O -i {resp1} -o {respout1} -e {esp} -t {qout1}".format(
              resp1=resp1_filename, respout1=respout1_filename,
              esp=esp_filename, qout1=qout1_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)

        #Stage 2
        qout2_filename = self.molname + '.qout_stage2'
        respout2_filename = self.molname + '.respout2'
        print("Using resp to perform stage 2 RESP fitting to generate: " + qout2_filename)

        cmd = "resp -O -i {resp2} -o {respout2} -e {esp} -q {qout1} -t {qout2}".format(
              resp2=resp2_filename, respout2=respout2_filename, esp=esp_filename,
              qout1=qout1_filename, qout2=qout2_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)
        
        # Write molecule with updated charges to mol2 file.
        mol2_filename = self.molname + ".mol2"
        print("Writing out the mol2 file with resp charge: " + mol2_filename)
        cmd = "antechamber -i {ac} -fi ac -o {mol2} -fo mol2 -c rc -cf {qout2}".format(
               ac=ac_filename, mol2=mol2_filename, qout2=qout2_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)

    def run(self):    
        # Run GAMESS to calculate ESP potential
        if not os.path.isdir(self.resp_scr_dir):
            print("Creating the scratch folder for RESP fitting: ", self.resp_scr_dir)
            os.mkdir(self.resp_scr_dir)

        print("Copying over the pdb file", self.pdbfile, " to ", self.resp_scr_dir)
        shutil.copy(os.path.join(self.rundir,self.pdbfile), self.resp_scr_dir)

        print("Navigating to the scratch folder for RESP fitting: ", self.resp_scr_dir)
        os.chdir(self.resp_scr_dir)

        self.runGAMESS()
        self.respFit()

        print("Navigating back to the folder for AutoSolvate run: ", self.rundir)
        os.chdir(self.rundir)
        
        mol2_filename = self.molname + ".mol2"
        print("Copying the generated mol2 file ", mol2_filename, " to rundir", self.rundir)
        shutil.copy(os.path.join(self.resp_scr_dir, mol2_filename), self.rundir)
