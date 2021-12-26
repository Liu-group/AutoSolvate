from openbabel import pybel
import getopt, sys, os
from openbabel import openbabel as ob
import subprocess


amber_solv_dict = {'water': [' ','TIP3PBOX '],
                   'methanol': ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
                   'chloroform': ['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
                   'nma': ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}

class solventBoxBuilder():
    r"""
    Solvated molecule in specified solvent.
    
    Parameters
    ----------
    solvent : str, Optional, default: 'water'
        Currently implemented solvents are: 'water', 'methanol', 'chloroform', 'nma', 'acetonitrile'
    slu_netcharge: int, Optional, default 0
        Charge of solute, the solvent box will be neutralized with Cl- and Na+ ions
    cube_size: float, Optional, default: 54
        Size of MM solvent box
    charge_method: str, Optional, default: "resp"
        Use 'resp' (quantum mechanical calculation needed) or 'bcc' to estimate partial charges
    slu_spinmult: int, Optional, default: 1
        Spinmultiplicity of solute
    outputFile: str, Optional, default='water_solvated'
        Filename-prefix for outputfiles
    srun_use: bool, Optional, default='False
        Run all commands with a srun prefix
    Returns
    -------
    None
        To run solvation, call build function.

    
    """


    def __init__(self, xyzfile, solvent='water', slu_netcharge=0, cube_size=54, 
            charge_method="resp", slu_spinmult=1, outputFile='water_solvated', 
            srun_use=False,gaussianexe=None, gaussiandir=None, amberhome=None, tolerance=2):
        self.xyz = xyzfile
        self.solute = pybel.readfile('xyz', xyzfile).__next__()
        self.slu_netcharge = slu_netcharge
        self.slu_spinmult = slu_spinmult
        # currently hard coded. Can be changed later to be determined automatically based on the density of the solute
        self.solvent = solvent
        self.tolerance=tolerance
        self.waterbox_size = 8.0
        # following are for custom organic solvent
        self.cube_size = cube_size # in angstrom
        self.slu_pos = self.cube_size/2.0
        self.slv_count = 210*8 # 210
        self.pbcbox_size = self.cube_size+2
        self.outputFile=outputFile
        self.srun_use=srun_use
        self.charge_method=charge_method
        self.gaussian_dir = gaussianexe
        self.gaussian_exe = gaussiandir
        self.amberhome = amberhome
        self.inputCheck()

    def inputCheck(self):
        if self.slu_spinmult > 1:
            if self.charge_method != "resp":
                print("Error: solute spin multiplicity: ", self.slu_spinmult, " charge method: ", self.charge_method)
                print("Error: atomic charge fitting for open-shell system only works for resp charge method")
                print("Error: exiting...")
                exit(1)
        if self.charge_method == "resp":
            if self.gaussian_exe == None:
                print("WARNING: Gaussian executable name is not specified for RESP charge fitting!")
                print("WARNING: Using g16 by default. If failed later, please rerun with the option -e specified!")
                self.gaussian_exe = 'g16'
            if self.gaussian_dir == None:
                print("WARNING: Gaussian executable directory is not specified for RESP charge fitting!")
                print("WARNING: Setting to default path: /opt/packages/gaussian/g16RevC.01/g16/")
                print("WARNING: If failed later, please rerun with the option -d specified!")
                self.gaussian_dir = '/opt/packages/gaussian/g16RevC.01/g16/'
                gspath = os.path.join(self.gaussian_dir, self.gaussian_exe)
                if not os.path.exists(gspath):
                    print("Error! Gaussian executable path",gspath)
                    print("Does not exist! Exiting...")
                    exit()
        if self.amberhome == None:
            print("WARNING: Amber home directory is not specified in input options")
            print("WARNING: Checking AMBERHOME environment virable...")
            cmd = ["echo", "$AMBERHOME"]
            print(cmd)
            proc=subprocess.Popen(cmd,
                           universal_newlines=True,
                           stdout=subprocess.PIPE)
            output=proc.communicate()[0]
            if len(output)<2:
                print("ERROR: AMBERHOME not defined")
                print("ERROR: Please set up AMBERHOME environment or specify through command line option -a")
                print("ERROR: Exiting...")
            else:
                print("WARNING: AMBERHOME detected: ", output)
                self.amberhome = output
        else:
            print("AMBERHOME path provided from input: ", self.amberhome)
            print("Validating path...")
            if not os.path.isdir(self.amberhome):
                print("ERROR: provided AMBERHOME path does not exist! Exiting...")

            print("Exporting AMBERHOME environment variable:")
            cmd = "export AMBERHOME=" + self.amberhome
            print(cmd)
            subprocess.call(cmd, shell=True)
            print("AMBERHOME environment variable export finished.")


    def getSolutePDB(self):
        r"""
        Convert xyz to pdb
        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        print("Converting xyz to pdb")
#        cmd = "babel " + self.xyz + " solute.pdb"
#        subprocess.call(cmd, shell=True)
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats("xyz", "pdb")
        obmol = self.solute.OBMol
        obConversion.WriteFile(obmol,'solute.xyz.pdb')
        # Change the residue name from the default UNL to SLU
        pdb1 = open('solute.xyz.pdb').readlines()
        pdb2 = open('solute.xyz.pdb','w')
        for line in pdb1:
            newline = line.replace('UNL','SLU')
            pdb2.write(newline)
        pdb2.close()

    def removeConectFromPDB(self):
        print("cleaning up solute.xyz.pdb")
        pdb1 = open('solute.xyz.pdb').readlines()
        pdb2 = open('solute.xyz.pdb','w')
        for line in pdb1:
            if 'CONECT' not in line:
                pdb2.write(line)
        pdb2.close()


            
    def getFrcmod(self):
        r"""
        Get partial charges and create frcmod

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Generate frcmod file for the solute")
        if self.charge_method == "resp":
            print("First generate the gaussian input file for RESP charge fitting")
            cmd1 ="$AMBERHOME/bin/antechamber -i solute.xyz.pdb -fi pdb -o gcrt.com -fo gcrt -gv 1 -ge solute.gesp  -s 2 -nc "+str(self.slu_netcharge) + " -m " + str(self.slu_spinmult)
            if self.srun_use:
                cmd1='srun -n 1 '+cmd1
            print(cmd1)
            subprocess.call(cmd1, shell=True)
            print("Then run Gaussian...")
            basedir=os.getcwd()
            if not os.path.isdir('tmp_gaussian'):
                os.mkdir('tmp_gaussian')

            cmd21="export GAUSS_EXEDIR=" + self.gaussian_dir + "; export GAUSS_SCRDIR="+basedir+"/tmp_gaussian; "
            #$PROJECT/TMP_GAUSSIAN;" #/expanse/lustre/projects/mit181/eh22/TMP_GAUSSIAN/;" #/scratch/$USER/$SLURM_JOBID ;"
            cmd22 = self.gaussian_exe + " < gcrt.com > gcrt.out"
            if self.srun_use:
                cmd22='srun -n 1 '+cmd22
            cmd2=cmd21+cmd22
            print(cmd2)
            subprocess.call(cmd2, shell=True)
            if not os.path.isfile('solute.gesp'):
                print("gaussian failed to generate solute.gesp")
                sys.stdout.flush()
                sys.exit()
            print("Gaussian ESP calculation done")
        if self.charge_method == "bcc":
           self.removeConectFromPDB()
        print("Then write out mol2")
        if self.charge_method == "resp":
            cmd3="$AMBERHOME/bin/antechamber -i solute.gesp -fi gesp -o solute.mol2 -fo mol2 -c resp -eq 2 -rn SLU"
        elif self.charge_method == "bcc":
            cmd3="$AMBERHOME/bin/antechamber -i solute.xyz.pdb -fi pdb -o solute.mol2 -fo mol2 -c bcc -eq 2 -rn SLU"
        if self.srun_use:
                 cmd3='srun -n 1 '+cmd3
        subprocess.call(cmd3, shell=True)
        print("Finally generate frcmod with parmchk2")
        cmd4 ="$AMBERHOME/bin/parmchk2 -i solute.mol2 -f mol2 -o solute.frcmod"
        if self.srun_use:
                cmd4='srun -n 1 '+cmd4
        subprocess.call(cmd4, shell=True)

    def getHeadTail(self):
        r"""
        Detect start and end of coordinates

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        pdb = open('solute.mol2').readlines()
        start =0
        end = 0
        for i in range(0,len(pdb)):
            if "@<TRIPOS>ATOM" in pdb[i]:
                start = i+1
            if "@<TRIPOS>BOND" in pdb[i]:
                end = i
        for i in range(start, end):
            atom =  pdb[i].split()[1]
            if "H" not in atom:
                self.head = atom
                break
        for i in reversed(range(start,end)):
            atom =  pdb[i].split()[1]
            if "H" not in atom:
                self.tail = atom
                break

    def writeTleapcmd1(self):
        r"""
        Create tleap input file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        self.getHeadTail()
        f = open("leap.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("loadamberparams solute.frcmod\n")
        f.write("SLU = loadmol2 solute.mol2\n")
        f.write("check SLU\n")
        f.write("set SLU head SLU.1."+self.head + "\n")
        f.write("set SLU tail SLU.1."+self.tail + "\n")
        f.write("saveoff SLU solute.lib\n")
        f.write("savepdb SLU solute.pdb\n")
        f.write("quit\n")
        f.close()
    
    def createLib(self):
        r"""
        Run tleap

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Now create the solute library file")
        self.writeTleapcmd1()
        cmd ="tleap -s -f leap.cmd > leap_savelib.log"
        if self.srun_use:
            cmd='srun -n 1 '+cmd
        subprocess.call(cmd, shell=True)
        if not os.path.isfile('solute.pdb'):
            print("gaussian failed to generate solute.pdb") 
            sys.stdout.flush()
            sys.exit()


    def writeTleapcmd_add_solvent(self):
        r"""
        Write tleap input file to add solvent

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if self.solvent in amber_solv_dict:
            print("Now add pre-equlibrated solvent box to the solute")
            self.getHeadTail()
            f = open("leap_add_solventbox.cmd","w")
            f.write("source leaprc.protein.ff14SB\n")
            f.write("source leaprc.gaff\n")
            f.write("source leaprc.water.tip3p\n")
            f.write(str(amber_solv_dict[str(self.solvent)][0]))        
            f.write("loadamberparams solute.frcmod\n")
            f.write("mol=loadmol2 solute.mol2\n")
            f.write("check mol\n")
            f.write("solvatebox mol " + (str(amber_solv_dict[str(self.solvent)][1])) + str(self.slu_pos) + " iso 0.8  #Solvate the complex with a cubic water box\n") 
            # Notice that we want to add the ion after solvation because we don't want the counter ion to be too close to solute
            if self.slu_netcharge != 0:
                if self.slu_netcharge > 0:
                    ion = 'Cl-'
                else:
                    ion = 'Na+'
                f.write("addIons2 mol " + ion + " 0\n")
                f.write("check mol\n")
            f.write("check mol\n")
            f.write("savepdb mol " + str(self.outputFile) + ".pdb\n")
            f.write("saveamberparm mol " + str(self.outputFile) + ".prmtop " + str(self.outputFile) + ".inpcrd\n")
            f.write("quit\n")
            f.close()


    def processPackmolPDB(self):
        r"""
        Convert file for CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        data=open('ch3cn_solvated.packmol.pdb').readlines()
        output=open('ch3cn_solvated.processed.pdb','w')
        this_resid = 1
        last_resid = 1
        for line in data:
            if 'ATOM' in line:
                last_resid = int(this_resid)
                this_resid = int(line[22:26])
            if last_resid != this_resid:
                output.write("TER\n")
            output.write(line)
        output.close()
    
    def packSLUSLV(self):
        r"""
        Write packmol input file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Now use packmol to pack the solute and solvent into a box")
        if self.solvent != 'acetonitrile':
            print("The type of solvent is not implemented yet")
            return
        else:
            packmol_inp = open('packmol.inp','w')
            packmol_inp.write("# All atoms from diferent molecules will be at least 2.0 Angstroms apart\n")
            packmol_inp.write("tolerance %s\n" % self.tolerance)
            packmol_inp.write("\n")
            packmol_inp.write("filetype pdb\n")
            packmol_inp.write("\n")
            packmol_inp.write("output ch3cn_solvated.packmol.pdb\n")
            packmol_inp.write("\n")
            packmol_inp.write("# add the solute\n")
            packmol_inp.write("structure solute.pdb\n")
            packmol_inp.write("   number 1\n")
            packmol_inp.write("   fixed " + " " + str(self.slu_pos) + " "+ str(self.slu_pos) + " " + str(self.slu_pos) + " 0. 0. 0.\n")
            packmol_inp.write("   centerofmass\n")
            packmol_inp.write("   resnumbers 2 \n")
            packmol_inp.write("end structure\n")
            packmol_inp.write("\n")
            packmol_inp.write("# add first type of solvent molecules\n")
            packmol_inp.write("structure ch3cn.pdb\n")
            packmol_inp.write("  number " + str(self.slv_count) + " \n")
            packmol_inp.write("  inside cube 0. 0. 0. " + str(self.cube_size) + " \n")
            packmol_inp.write("  resnumbers 2 \n")
            packmol_inp.write("end structure\n")
            packmol_inp.close()

            cmd ="packmol < packmol.inp > packmol.log"
            subprocess.call(cmd, shell=True)
            if self.srun_use:
                            cmd='srun -n 1 '+cmd
            self.processPackmolPDB()

    def writeTleapcmd_ch3cn_solvated(self):
        r"""
        Write tleap file for CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        f = open("leap_packmol_solvated.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
        f.write("loadamberparams ch3cn.frcmod\n")
        f.write("loadAmberPrep ch3cn.prep\n")
        f.write("loadamberparams solute.frcmod\n")
        f.write("loadoff solute.lib\n")
        f.write("\n")
        f.write("SYS = loadpdb ch3cn_solvated.processed.pdb\n")
        f.write("check SYS\n")
        f.write("\n")
        if self.slu_netcharge > 0:
            ion = 'Cl-'
        else:
            ion = 'Na+'
        if self.slu_netcharge != 0:
            if self.slu_netcharge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            f.write("addIons2 SYS " + ion + " 0\n")
            f.write("check SYS\n")
        f.write("# set the dimension of the periodic box\n")
        f.write("set SYS box {" + str(self.pbcbox_size) +", " + str(self.pbcbox_size) + ", " + str(self.pbcbox_size) + "}\n")
        f.write("\n")
        f.write("saveamberparm SYS "+str(self.outputFile)+".prmtop "+str(self.outputFile)+".inpcrd #Save AMBER topology and coordinate files\n")
        f.write("savepdb SYS "+str(self.outputFile)+".pdb\n")
        f.write("quit\n")
        f.close()

    def createAmberParm(self):
        r"""
        Generate Amber parameters with tleap

        Parameters
        ---------
        None

        Returns
        ---------
        None
            Creates files in current directory.

        """
        print("Generate Amber parameters for the solvated system")
        if self.solvent in amber_solv_dict:
            self.writeTleapcmd_add_solvent()
            cmd ="tleap -s -f leap_add_solventbox.cmd > leap_add_solventbox.log"
            if self.srun_use:
                            cmd='srun -n 1 '+cmd
            subprocess.call(cmd, shell=True)
        elif self.solvent == 'acetonitrile':
            self.packSLUSLV()
            self.writeTleapcmd_ch3cn_solvated()
            cmd ="tleap -s -f leap_packmol_solvated.cmd > leap_packmol_solvated.log"
            if self.srun_use:
                            cmd='srun -n 1 '+cmd
            subprocess.call(cmd, shell=True)
        else:
            print('solvent not supported')

    def build(self):
        r"""
        Run solventBoxBuilder
         
        Parameters
        ---------
        None
       
        Returns
        ---------
        None
	    Creates files in current directory.

        """
        self.getSolutePDB()
        self.getFrcmod()
        self.createLib()
        self.createAmberParm()
        print("The script has finished successfully")

def startboxgen(argumentList):
    r"""
    Wrap function that parses commandline options,
    adds solvent box to a given solute,
    and generates related force field parameters.
    
    Parameters
    ----------
    argumentList: list
       The sist contains the command line options to specify solute, solvent, and other options
       related structure and force field parameter generation.

       Flag definitions:

       -m, --main, solute xyz file
       -s, --solvent, name of solvent (water, methanol, chloroform, nma)
       -o, --output, prefix of the output file names
       -c, --charge, formal charge of solute 
       -u, --spinmultiplicity, spin multiplicity of solute 
       -g, --chargemethod, name of charge fitting method (bcc, resp)
       -b, --cubesize, size of solvent cube in angstroms 
       -r, --srunuse, option to use slurm to submit as a job
       -e, --gaussianexe, name of the Gaussian quantum chemistry package executable used to generate electrostatic potential needed for RESP charge fitting
       -d, --gaussiandir, path to the Gaussian package.
       -a, --amberhome, path to the AMBER molecular dynamics package root directory. Definition of the environment variable $AMBERHOME.


    Returns
    -------
    None
    """
    print(argumentList)
    options = "m:s:o:c:k:b:g:u:r:e:d:a:t"
    long_options = ["main", "solvent", "output", "charge", "cubesize", "chargemethod", "spinmultiplicity", "srunuse","gaussianexe", "gaussiandir", "amberhome", "tolerance"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    solutexyz=""
    solvent='water'
    slu_netcharge=0
    cube_size=54
    charge_method="bcc"
    slu_spinmult=1
    outputFile='water_solvated'
    srun_use=False
    amberhome=None
    gaussianexe=None
    gaussiandir=None
    tolerance=2
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-m", "-main"):
            print ("Main/solutexyz", currentValue)
            solutexyz=str(currentValue)     
        elif currentArgument in ("-s", "-solvent"):
            print ("Solvent:", currentValue)
            solvent=str(currentValue)
        elif currentArgument in ("-o", "-output"):
            print ("Output:", currentValue)
            outputFile=str(currentValue)
        elif currentArgument in ("-c", "-charge"):
            print ("Charge:", currentValue)
            slu_netcharge=int(currentValue)
        elif currentArgument in ("-b", "-cubesize"):
            print ("Cubesize:", currentValue)
            cube_size=float(currentValue)
        elif currentArgument in ("-g", "-chargemethod"):
            print ("Chargemethod:", currentValue)
            charge_method=str(currentValue)
        elif currentArgument in ("-u", "-spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            slu_spinmult=int(currentValue)
        elif currentArgument in ("-r", "-srunuse"):
            print("usign srun")
            srun_use=True
        elif currentArgument in ("-e","gaussianexe"):
            print("Gaussian executable name:", currentValue)
            gaussianexe = currentValue
        elif currentArgument in ("-d","gaussiandir"):
            print("Gaussian package directory:", currentValue)
            gaussiandir = currentValue
        elif currentArgument in ("-a","amberhome"):
            print("Amber home directory:", currentValue)
            amberhome = currentValue
        elif currentArgument in ("-t", "tolerance"):
            print("Tolerance for Packmol", currentValue)
            tolerance = currentValue

    if solutexyz == "":
        print("Error! Solute xyzfile must be provided!\nExiting...")
        exit()
    elif not os.path.exists(solutexyz):
        print("Error! Solute xyzfile path ",solutexyz, " does not exist!\nExiting...")
        exit()

    try:
        pybel.readfile('xyz', solutexyz).__next__()
    except:
        print("Error! Solute xyzfile format issue!")
        print(solutexyz," cannot be opened with openbabel.\n Exiting...")
        exit()
     
    builder = solventBoxBuilder(solutexyz, solvent, slu_netcharge, cube_size, charge_method, 
                                slu_spinmult, outputFile, srun_use=srun_use, 
                                gaussianexe=gaussianexe, gaussiandir=gaussiandir, amberhome=amberhome, tolerance=tolerance)
    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startboxgen(argumentList)
