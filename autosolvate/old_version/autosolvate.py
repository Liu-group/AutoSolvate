from openbabel import pybel
import getopt, sys, os
from openbabel import openbabel as ob
import subprocess
import pkg_resources
import numpy as np

from .terachem_docker import *


amber_solv_dict = {'water':     [' ','TIP3PBOX '],
                   'methanol':  ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
                   'chloroform':['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
                   'nma':       ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}

custom_solv_dict = {'acetonitrile':'ch3cn'}

class AmberParamsBuilder():
    """ 
    This class handles the Amber parameter creation for one single molecule.
    1. Generate standard pdb
    2. AnteChamber or Gaussian charge fitting
    3. Tleap create Lib
    """
    def __init__(self, xyzfile:str, name = "", resname = "", charge = 0, spinmult = 1, charge_method="resp", 
        outputFile="", srun_use=False, gaussianexe=None, gaussiandir=None, amberhome=None, use_terachem = False):

        self.xyz = xyzfile
        self.xyzformat = os.path.splitext(self.xyz)[-1][1:]
        print(f"{self.xyzformat.upper()} file detected by checking the extension.")
        if not name:
            self.name = os.path.basename(xyzfile)
            self.name = os.path.splitext(self.xyz)[0]
        else:
            self.name = name

        self.resname = self.name[:3].upper() if not resname else resname[:3].upper()
        print("Set residue name of component " + self.name + " as " + self.resname)

        self.molecule = pybel.readfile(self.xyzformat, self.xyz).__next__()
        # following are for custom organic solvent
        self.outputFile = outputFile
        self.srun_use = srun_use
        self.charge = charge
        self.spinmult = spinmult
        self.charge_method = charge_method
        self.gaussian_dir = gaussiandir
        self.gaussian_exe = gaussianexe
        self.amberhome = amberhome
        self.use_terachem = use_terachem
        self.inputCheck()

    def inputCheck(self):
        if not os.path.exists(self.xyz):
            print("Error: structure file " + self.xyz + "not found")
            exit(1)
        if self.spinmult > 1:
            if self.charge_method != "resp":
                print("Error: spin multiplicity: ", self.spinmult, " charge method: ", self.charge_method)
                print("Error: atomic charge fitting for open-shell system only works for resp charge method")
                print("Error: exiting...")
                exit(1)
        if self.charge_method == "resp" and self.use_terachem:
            check_terachem()
        elif self.charge_method == "resp" and not self.use_terachem:
            if self.gaussian_exe == None:
                print("WARNING: Gaussian executable name is not specified for RESP charge fitting!")
                print("WARNING: Using g16 by default. If failed later, please rerun with the option -e specified!")
                self.gaussian_exe = 'g16'
            if self.gaussian_dir == None or not os.path.exists(self.gaussian_dir):
                print("WARNING: Gaussian executable directory is not specified for RESP charge fitting!")
                print("WARNING: Setting to default path: /opt/packages/gaussian/g16RevC.01/g16/")
                print("WARNING: If failed later, please rerun with the option -d specified!")
                self.gaussian_dir = '/opt/packages/gaussian/g16RevC.01/g16/'
                gspath = os.path.join(self.gaussian_dir, self.gaussian_exe)
                if not os.path.exists(gspath):
                    print("Error! Gaussian executable path",gspath)
                    print("Does not exist! Exiting....")
                    exit()
        if self.amberhome == None:
            print("WARNING: Amber home directory is not specified in input options")
            print("WARNING: Checking AMBERHOME environment variable...")
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

    def getPDB(self):
        r"""
        Convert input structure to pdb. 
        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        print("Converting to pdb")
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats(self.xyzformat, "pdb")
        obmol = self.molecule.OBMol
        obConversion.WriteFile(obmol, f"{self.name}.pdb")
        # Change the residue name from the default UNL to SLV
        res = obmol.GetResidue(0)
        res.SetName(self.resname)
        obConversion.WriteFile(obmol, f"{self.name}.pdb")

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
        pdb = open(f'{self.name}.mol2').readlines()
        start = 0
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

    def removeConectFromPDB(self):
        print(f"cleaning up {self.name}.pdb")
        pdb1 = open(f'{self.name}.pdb').readlines()
        pdb2 = open(f'{self.name}.pdb','w')
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
        print(f"Generate frcmod file for {self.name}")
        if self.charge_method == "resp":
            if self.use_terachem:
                resp_terachem(self.name + ".pdb", self.resname, self.charge, self.spinmult)
                cmd4 = f"$AMBERHOME/bin/parmchk2 -i {self.name}.mol2 -f mol2 -o {self.name}.frcmod"
                if self.srun_use:
                    cmd4 = 'srun -n 1 ' + cmd4
                subprocess.call(cmd4, shell=True)
                return 
            print("First generate the gaussian input file for RESP charge fitting")
            cmd1 = f"$AMBERHOME/bin/antechamber -i {self.name}.pdb -fi pdb -o gcrt.com -fo gcrt -gv 1 -ge {self.name}.gesp  -s 2 -nc {self.charge} -m {self.spinmult}"
            if self.srun_use:
                cmd1 = 'srun -n 1 ' + cmd1
            print(cmd1)
            subprocess.call(cmd1, shell=True)
            print("Then run Gaussian...")
            basedir = os.getcwd()
            if not os.path.isdir('tmp_gaussian'):
                os.mkdir('tmp_gaussian')

            cmd21 = "export GAUSS_EXEDIR=" + self.gaussian_dir + "; export GAUSS_SCRDIR="+basedir+"/tmp_gaussian; "
            #$PROJECT/TMP_GAUSSIAN;" #/expanse/lustre/projects/mit181/eh22/TMP_GAUSSIAN/;" #/scratch/$USER/$SLURM_JOBID ;"
            cmd22 = self.gaussian_exe + " < gcrt.com > gcrt.out"
            if self.srun_use:
                cmd22 = 'srun -n 1 ' + cmd22
            cmd2 = cmd21 + cmd22
            print(cmd2)
            subprocess.call(cmd2, shell=True)
            if not os.path.isfile(f'{self.name}.gesp'):
                print(f"gaussian failed to generate {self.name}.gesp")
                raise Exception
                sys.stdout.flush()
                sys.exit()
            print("Gaussian ESP calculation done.")
        print("Do charge generation.")
        pl = -1
        if len(self.molecule.atoms) >= 100:
            pl = 15
        if self.charge_method == "bcc":
            self.removeConectFromPDB()
        if self.charge_method == "resp":
            cmd3 = f"$AMBERHOME/bin/antechamber -i {self.name}.gesp -fi gesp -o {self.name}.mol2 -fo mol2 -c resp -eq 2 -rn {self.resname} -pl {pl} -nc {self.charge} -m {self.spinmult}"
        elif self.charge_method == "bcc":
            cmd3 = f"$AMBERHOME/bin/antechamber -i {self.name}.pdb -fi pdb -o {self.name}.mol2 -fo mol2 -c bcc -eq 2 -rn {self.resname} -pl {pl} -nc {self.charge} -m {self.spinmult}"
        if self.srun_use:
            cmd3 = 'srun -n 1 ' + cmd3
        subprocess.call(cmd3, shell=True)
        print("Finally generate frcmod with parmchk2")
        cmd4 = f"$AMBERHOME/bin/parmchk2 -i {self.name}.mol2 -f mol2 -o {self.name}.frcmod"
        if self.srun_use:
            cmd4 = 'srun -n 1 ' + cmd4
        subprocess.call(cmd4, shell=True)

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
        f = open(f"leap_{self.name}.cmd","w")
        f.write(f"source leaprc.protein.ff14SB\n")
        f.write(f"source leaprc.gaff\n")
        f.write(f"loadamberparams {self.name}.frcmod\n")
        f.write(f"{self.resname} = loadmol2 {self.name}.mol2\n")
        f.write(f"check {self.resname}\n")
        f.write(f"set {self.resname} head {self.resname}.1."+self.head + "\n")
        f.write(f"set {self.resname} tail {self.resname}.1."+self.tail + "\n")
        f.write(f"saveoff {self.resname} {self.name}.lib\n")
        f.write(f"savepdb {self.resname} {self.name}.pdb\n")
        f.write(f"saveamberparm {self.resname} {self.name}.prmtop {self.name}.inpcrd\n")
        f.write("quit\n")
        f.close()

    def convert2charmm(self):
        subprocess.run(f"amb2chm_psf_crd.py -p {self.name}.prmtop -c {self.name}.inpcrd -f {self.name}.psf -d {self.name}.crd -b {self.name}.crd.pdb", shell = True)
        amberhome = self.amberhome
        if amberhome.find("AMBERHOME") != -1:
            amberhome = os.getenv("AMBERHOME")
        f = open(f"amb2chm_{self.name}.in", "w")
        f.write(f"{amberhome}/dat/leap/parm/gaff.dat\n")
        # f.write(f"{amberhome}/dat/leap/parm/frcmod.ff19SB\n")
        # amb2chm_par.py parmed.exceptions.ParameterError: Could not understand nonbond parameter line [CMAP]
        f.write(f"{self.name}.frcmod")
        f.close()
        subprocess.run(f"amb2chm_par.py -i amb2chm_{self.name}.in -o {self.name}.prm", shell = True)
    
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
        print(f"Now create the library file for {self.name}")
        self.writeTleapcmd1()
        cmd = f"tleap -s -f leap_{self.name}.cmd > leap_{self.name}_savelib.log"
        if self.srun_use:
            cmd='srun -n 1 ' + cmd
        subprocess.call(cmd, shell=True)
        if not os.path.isfile(f'{self.name}.pdb'):
            print(f"gaussian failed to generate {self.name}.pdb") 
            sys.stdout.flush()
            sys.exit()

    def build(self):
        r"""
        Run ParameterBuilder
         
        Parameters
        ---------
        None
       
        Returns
        ---------
        None
	    Creates files in current directory.

        """         
        self.getPDB()
        self.getFrcmod()
        self.createLib()
        print("The script has finished successfully")

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

    def __init__(self, xyzfile:str, solvent='water', slu_netcharge=0, cube_size=54, 
            charge_method="resp", slu_spinmult=1, outputFile="", 
            srun_use=False, gaussianexe=None, gaussiandir=None, amberhome=None, use_terachem = False,
            closeness=0.8, solvent_off="", solvent_frcmod="",
            slu_count=1, slv_count=210*8, slv_generate=False, slv_xyz=""):
        self.xyz = xyzfile
        self.xyzformat = os.path.splitext(self.xyz)[-1][1:]
        print(f"{self.xyzformat.upper()} file detected by checking the extension.")
        self.solute = pybel.readfile(self.xyzformat, self.xyz).__next__()
        # currently hard coded. Can be changed later to be determined automatically based on the density of the solute
        self.solvent = solvent
        if closeness=='automated':
          if self.solvent=='acetonitrile':
            self.closeness=1.88
          elif self.solvent=='water':
            self.closeness=0.50
          elif self.solvent=='methanol':
            self.closeness=0.60
          elif self.solvent=='nma':
            self.closeness=0.58
          elif self.solvent=='chloroform':
            self.closeness=0.58
          else: 
            print("Warning: solvent not supported for automated closeness determination")           
        else:
          self.closeness = closeness
        self.waterbox_size = 8.0
        # following are for custom organic solvent
        self.cube_size = cube_size # in angstrom
        self.pbcbox_size = self.cube_size+2
        self.outputFile=outputFile
        if self.outputFile == "":
            self.outputFile = self.solvent + "_solvated"
        self.srun_use=srun_use
        self.charge_method=charge_method
        self.gaussian_dir = gaussiandir
        self.gaussian_exe = gaussianexe
        self.amberhome = amberhome
        self.use_terachem = use_terachem
        self.solvent_off=solvent_off
        self.solvent_frcmod=solvent_frcmod
        self.is_custom_solvent = False

        self.slu_name = os.path.splitext(os.path.basename(xyzfile))[0]
        self.slu_netcharge = slu_netcharge
        self.slu_spinmult = slu_spinmult
        self.slu_count = slu_count
        self.slu_pos = cube_size/2.0

        self.slv_name = self.solvent
        self.slv_count = slv_count # 210
        self.slv_generate = slv_generate
        self.slv_xyz = slv_xyz
        self.slv_xyzformat = os.path.splitext(self.slv_xyz)[-1]
        self.inputCheck()

    def inputCheck(self):
        if self.slu_spinmult > 1:
            if self.charge_method != "resp":
                print("Error: solute spin multiplicity: ", self.slu_spinmult, " charge method: ", self.charge_method)
                print("Error: atomic charge fitting for open-shell system only works for resp charge method")
                print("Error: exiting...")
                exit(1)
        if self.charge_method == "resp" and self.use_terachem:
            check_terachem()
        elif self.charge_method == "resp" and not self.use_terachem:
            assert self.use_terachem
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
                    print("Does not exist! Exiting....")
                    exit()
        if self.amberhome == None:
            print("WARNING: Amber home directory is not specified in input options")
            print("WARNING: Checking AMBERHOME environment variable...")
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
        # check whether requested custom solvent is available
        if self.solvent not in amber_solv_dict.keys() and self.solvent not in custom_solv_dict.keys():
            print("Requested solvent name not contained in AutoSolvate.")
            print("Checking available of custom solvent .frcmod and .off files")
            if not self.slv_generate:
                if len(self.solvent_frcmod) == 0:
                    print("ERROR! Custom solvent .frcmod file is not provided!")
                    exit()
                elif len(self.solvent_off) == 0:
                    print("ERROR! Custom solvent .off library file is not provided!")
                    exit()
                elif not os.path.exists(self.solvent_frcmod):
                    print("ERROR! Custom solvent .frcmod file ", self.solvent_frcmod,
                        "ERROR! does not exist!")
                    exit()
                elif not os.path.exists(self.solvent_off):
                    print("ERROR! Custom solvent .off library file ", self.solvent_frcmod,
                        "ERROR! does not exist!")
                    exit()
                else:
                    self.is_custom_solvent = True
            elif self.slv_generate:
                if len(self.solvent_frcmod) == 0:
                    print("Custom solvent .frcmod file is not provided! Will generate from GAFF instead.")
                elif len(self.solvent_off) == 0:
                    print("Custom solvent .off library file  is not provided! The box will generate from packmol instead.")
                elif not os.path.exists(self.solvent_frcmod):
                    print("Custom solvent .frcmod file ", self.solvent_frcmod,
                        "does not exist! Will generate from GAFF instead.")
                elif not os.path.exists(self.solvent_off):
                    print("Custom solvent .off library file ", self.solvent_frcmod,
                        "does not exist! The box will generate from packmol instead.")
                else:
                    print("Solvent parameter generation is disabled as .frcmod and .off are provided.")
                    self.is_custom_solvent = True
                    self.slv_generate = False
            if self.slv_generate:
                if not self.slv_xyz or not os.path.exists(self.slv_xyz):
                    print("ERROR! The coordinate file for custom solvent must be provided when generate parameter from forcefield.")
                    exit()
                cmd = "which packmol"
                proc=subprocess.Popen(cmd, shell = True,
                           universal_newlines=True,
                           stdout=subprocess.PIPE)
                output=proc.communicate()[0]
                if not output:
                    print("ERROR! Cannot find packmol in $PATH. packmol must be installed when using custom solvent and multiple solute.")
                    exit()
            
        else:
            print("Solvent parameter generation is disabled as requested solvent is contained in Autosolvate.")
            self.slv_generate = False
        

        if len(self.solute.atoms) > 100:
            print("The solute molecule has over 100 atoms. Will adjust -pl to 15")
         
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
            f = open("leap_add_solventbox.cmd","w")
            f.write("source leaprc.protein.ff14SB\n")
            f.write("source leaprc.gaff\n")
            f.write("source leaprc.water.tip3p\n")
            f.write(str(amber_solv_dict[str(self.solvent)][0]))
            if os.path.exists(f"{self.slu_name}.frcmod"):
                f.write(f"loadamberparams {self.slu_name}.frcmod\n")
            f.write(f"loadoff {self.slu_name}.lib\n")
            f.write(f"mol=loadmol2 {self.slu_name}.mol2\n")
            f.write("check mol\n")
            f.write("solvatebox mol " + (str(amber_solv_dict[str(self.solvent)][1])) + str(self.slu_pos) + " iso "+str(self.closeness)+"  #Solvate the complex with a cubic solvent box\n") 
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

    def writeTleapcmd_add_solvent_custom(self):
        r"""
        Write tleap input file to add solvent from user provided .off and .frcmod files

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Now add custom pre-equlibrated solvent box to the solute")
        f = open("leap_add_solventbox.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n")
        f.write("loadoff " + self.solvent_off + "\n")
        f.write("loadamberparams " + self.solvent_frcmod + "\n")
        if os.path.exists(f"{self.slu_name}.frcmod"):
            f.write(f"loadamberparams {self.slu_name}.frcmod\n")
        f.write(f"mol=loadmol2 {self.slu_name}.mol2\n")
        f.write("check mol\n")
        f.write("solvatebox mol " + self.solvent + " " 
                + str(self.slu_pos) + " iso 0.8  #Solvate the complex with a cubic solvent box\n") 
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
        Convert file for custom solvents like CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        solvPrefix = self.solvent if self.slv_generate else custom_solv_dict[self.solvent]
        
        data=open(solvPrefix + '_solvated.packmol.pdb').readlines()
        output=open(solvPrefix + '_solvated.processed.pdb','w')
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
        Write packmol input file for custom solvents

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Now use packmol to pack the solute and solvent into a box")
        if not self.slv_generate:
            if self.solvent not in custom_solv_dict.keys():
                print("The type of solvent is not implemented yet")
                return
            else:
                # Solvent pdb file stored in our package data folder
                # Path is too long for Packmol to recognize. Copy to current rundir 
                solvPrefix = custom_solv_dict[self.solvent]
                solvent_pdb = solvPrefix+'.pdb'
                solvent_pdb_origin = pkg_resources.resource_filename('autosolvate', 
                        os.path.join('data/', solvPrefix, solvent_pdb))
                subprocess.call(['cp',solvent_pdb_origin,solvent_pdb])
        elif self.slv_generate:
            solvPrefix = self.solvent
            solvent_pdb = self.slv_name+".pdb"

        output_pdb = solvPrefix + "_solvated.packmol.pdb"

        packmol_inp = open('packmol.inp','w')
        packmol_inp.write("# All atoms from diferent molecules will be at least %s Angstroms apart\n" % self.closeness)
        packmol_inp.write("tolerance %s\n" % self.closeness)
        packmol_inp.write("\n")
        packmol_inp.write("filetype pdb\n")
        packmol_inp.write("\n")
        packmol_inp.write("output " + output_pdb + "\n")
        packmol_inp.write("\n")
        packmol_inp.write("# add the solute\n")
        if self.slu_count == 1:
            packmol_inp.write(f"structure {self.slu_name}.pdb\n")
            packmol_inp.write("   number 1\n")
            packmol_inp.write("   fixed " + " " + str(self.slu_pos) + " "+ str(self.slu_pos) + " " + str(self.slu_pos) + " 0. 0. 0.\n")
            packmol_inp.write("   centerofmass\n")
            packmol_inp.write("   resnumbers 2 \n")
            packmol_inp.write("end structure\n")
            packmol_inp.write("\n")
        elif self.slu_count > 1:
            packmol_inp.write("# add the solute\n")
            packmol_inp.write(f"structure {self.slu_name}.pdb\n")
            packmol_inp.write("    number %d\n" % self.slu_count)
            packmol_inp.write("    inside cube 0. 0. 0. " + str(self.cube_size) + " \n")
            packmol_inp.write("    resnumbers 2 \n")
            packmol_inp.write("end structure\n")
            packmol_inp.write("\n")

        packmol_inp.write("# add first type of solvent molecules\n")
        packmol_inp.write("structure "+ solvent_pdb + "\n")
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

    def writeTleapcmd_custom_solvated(self):
        r"""
        Write tleap file for custom solvents like CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if self.slv_generate:
            solvPrefix = self.solvent
            solvent_frcmod = self.slv_name+".frcmod"
            solvent_frcmod_path = os.path.join(os.getcwd(), solvent_frcmod)
            solvent_prep = self.slv_name+".prep"
            solvent_prep_path = os.path.join(os.getcwd(), solvent_prep)
            solvent_mol2 = self.slv_name+".mol2"
            solvent_mol2_path = os.path.join(os.getcwd(), solvent_mol2)
            solvent_lib = self.slv_name+".lib"
            solvent_lib_path = os.path.join(os.getcwd(), solvent_lib)
        else:
            solvPrefix = custom_solv_dict[self.solvent]
            solvent_frcmod = solvPrefix+'.frcmod'
            solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data',solvPrefix,solvent_frcmod))

            solvent_prep = solvPrefix+'.prep'
            solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data',solvPrefix,solvent_prep))
            solvent_lib_path = ""
            solvent_mol2_path = ""

        f = open("leap_packmol_solvated.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
        f.write("loadamberparams " + solvent_frcmod_path + "\n")
        for command, fpath in zip(["loadamberprep", "loadoff", "loadmol2"], [solvent_prep_path, solvent_lib_path, solvent_mol2_path]):
            print(command, fpath)
            if not fpath:
                continue
            if os.path.exists(fpath):
                f.write(command + " " + fpath + "\n")
                break
        else:
            print("ERROR: no solvent amberprep or library file!")
            exit()
        if os.path.exists(f"{self.slu_name}.frcmod"):
            f.write(f"loadamberparams {self.slu_name}.frcmod\n")
        f.write(f"loadoff {self.slu_name}.lib\n")
        f.write("\n")
        f.write("SYS = loadpdb " + solvPrefix + "_solvated.processed.pdb\n")
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
        elif self.solvent in custom_solv_dict.keys() or self.slv_generate:
            print("Custom solvent without solvent box, will use packmol to generate the solvent box.")
            self.packSLUSLV()
            self.writeTleapcmd_custom_solvated()
            cmd ="tleap -s -f leap_packmol_solvated.cmd > leap_packmol_solvated.log"
            if self.srun_use:
                cmd='srun -n 1 '+cmd
            subprocess.call(cmd, shell=True)
        elif self.is_custom_solvent:
            self.writeTleapcmd_add_solvent_custom()  
            cmd ="tleap -s -f leap_add_solventbox.cmd > leap_add_solventbox.log"
            if self.srun_use:
                cmd='srun -n 1 '+cmd
            subprocess.call(cmd, shell=True)
        else:
            print('solvent not supported')
            exit(1)

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
        solutebuilder = AmberParamsBuilder(
            xyzfile= self.xyz,
            name = self.slu_name,
            resname = "SLU",
            charge = self.slu_netcharge,
            spinmult= self.slu_spinmult,
            charge_method=self.charge_method,
            outputFile="solutegen.out",
            srun_use=self.srun_use,
            gaussianexe=self.gaussian_exe,
            gaussiandir=self.gaussian_dir,
            amberhome=self.amberhome,
            use_terachem=self.use_terachem
        )
        solutebuilder.build()

        if self.slv_generate:
            solventbuilder = AmberParamsBuilder(
                xyzfile = self.slv_xyz,
                name = self.slv_name,
                resname = "SLV",
                charge = 0,
                spinmult = 1,
                charge_method=self.charge_method,
                outputFile="solventgen.out",
                srun_use=self.srun_use,
                gaussianexe=self.gaussian_exe,
                gaussiandir=self.gaussian_dir, 
                amberhome=self.amberhome,
                use_terachem=self.use_terachem
            )
            solventbuilder.build()
            
        self.createAmberParm()
        print("The script has finished successfully")

def startboxgen(argumentList):
    r"""
    Wrap function that parses command line options for autosolvate boxgen,
    adds solvent box to a given solute,
    and generates related force field parameters.
    
    Parameters
    ----------
    argumentList: list
       The list contains the command line options to specify solute, solvent, and other options
       related to structure and force field parameter generation.

       Command line option definitions
         -m, --main  solute xyz file
         -s, --solvent  name of solvent (water, methanol, chloroform, nma)
         -o, --output  prefix of the output file names
         -c, --charge  formal charge of solute
         -u, --spinmultiplicity  spin multiplicity of solute
         -g, --chargemethod  name of charge fitting method (bcc, resp)
         -b, --cubesize  size of solvent cube in angstroms
         -r, --srunuse  option to run inside a slurm job
         -e, --gaussianexe  name of the Gaussian quantum chemistry package executable used to generate electrostatic potential needed for RESP charge fitting
         -d, --gaussiandir  path to the Gaussian package
         -a, --amberhome  path to the AMBER molecular dynamics package root directory. Definition of the environment variable $AMBERHOME
         -t, --closeness  Solute-solvent closeness setting, for acetonitrile tolerance parameter in packmol in Å, for water, methanol, nma, chloroform the scaling factor in tleap, setting to 'automated' will automatically set this parameter based on solvent.
         -l, --solventoff  path to the custom solvent .off library file. Required if the user want to use some custom solvent other than the 5 solvents contained in AutoSolvate (TIP3P water, methanol, NMA, chloroform, MeCN)
         -p, --solventfrcmod  path to the custom solvent .frcmod file. Required if the user wants to use some custom solvent other than the 5 solvents contained in AutoSolvate.
         -h, --help  short usage description

    Returns
    -------
    None
        Generates the structure files and save as ```.pdb```. Generates the MD parameter-topology and coordinates files and saves as ```.prmtop``` and ```.inpcrd```
    """
    #print(argumentList)
    options = "hm:s:o:c:b:g:u:re:d:a:t:l:p:"
    long_options = ["help", "main", "solvent", "output", "charge", "cubesize", "chargemethod", "spinmultiplicity", "srunuse","gaussianexe", "gaussiandir", "amberhome", "closeness","solventoff","solventfrcmod"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    solutexyz=""
    solvent='water'
    slu_netcharge=0
    cube_size=54
    charge_method="bcc"
    slu_spinmult=1
    outputFile=""
    srun_use=False
    amberhome=None
    gaussianexe=None
    gaussiandir=None
    closeness=0.8
    solvent_off=""
    solvent_frcmod=""
    #print(arguments)
    #print(values)
    for currentArgument, currentValue in arguments:
        if  currentArgument in ("-h", "--help"):
            print('Usage: autosolvate boxgen [OPTIONS]')
            print('  -m, --main                 solute xyz file')
            print('  -s, --solvent              name of solvent')
            print('  -o, --output               prefix of the output file names')
            print('  -c, --charge               formal charge of solute')
            print('  -u, --spinmultiplicity     spin multiplicity of solute')
            print('  -g, --chargemethod         name of charge fitting method (bcc, resp)')
            print('  -b, --cubesize             size of solvent cube in angstroms')
            print('  -r, --srunuse              option to run inside a slurm job')
            print('  -e, --gaussianexe          name of the Gaussian quantum chemistry package executable')
            print('  -d, --gaussiandir          path to the Gaussian package')
            print('  -a, --amberhome            path to the AMBER molecular dynamics package root directory')
            print('  -t, --closeness            Solute-solvent closeness setting')
            print('  -l, --solventoff           path to the custom solvent .off library file')
            print('  -p, --solventfrcmod        path to the custom solvent .frcmod file')
            print('  -h, --help                 short usage description')
            exit()
        elif currentArgument in ("-m", "--main"):
            print ("Main/solutexyz", currentValue)
            solutexyz=str(currentValue)     
        elif currentArgument in ("-s", "--solvent"):
            print ("Solvent:", currentValue)
            solvent=str(currentValue)
        elif currentArgument in ("-o", "--output"):
            print ("Output:", currentValue)
            outputFile=str(currentValue)
        elif currentArgument in ("-c", "--charge"):
            print ("Charge:", currentValue)
            slu_netcharge=int(currentValue)
        elif currentArgument in ("-b", "--cubesize"):
            print ("Cubesize:", currentValue)
            cube_size=float(currentValue)
        elif currentArgument in ("-g", "--chargemethod"):
            print ("Chargemethod:", currentValue)
            charge_method=str(currentValue)
        elif currentArgument in ("-u", "--spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            slu_spinmult=int(currentValue)
        elif currentArgument in ("-r", "--srunuse"):
            print("usign srun")
            srun_use=True
        elif currentArgument in ("-e","--gaussianexe"):
            print("Gaussian executable name:", currentValue)
            gaussianexe = currentValue
        elif currentArgument in ("-d","--gaussiandir"):
            print("Gaussian package directory:", currentValue)
            gaussiandir = currentValue
        elif currentArgument in ("-a","--amberhome"):
            print("Amber home directory:", currentValue)
            amberhome = currentValue
        elif currentArgument in ("-t", "--closeness"):
            print("Solute-Solvente closeness parameter", currentValue)
            closeness = currentValue
        elif currentArgument in ("-l", "--solventoff"):
            print("Custom solvent .off library path:", currentValue)
            solvent_off = currentValue
        elif currentArgument in ("-p", "--solventfrcmod"):
            print("Custom solvent .frcmmod file path:", currentValue)
            solvent_frcmod = currentValue

    if solutexyz == "":
        print("Error! Solute xyzfile must be provided!\nExiting...")
        exit()
    elif not os.path.exists(solutexyz):
        print("Error! Solute xyzfile path ", solutexyz, " does not exist!\nExiting...")
        exit()

    try:
        _, ext = os.path.splitext(solutexyz)
        pybel.readfile(ext[1:], solutexyz).__next__()
    except:
        print("Error! Solute structure file format issue!")
        print(solutexyz," cannot be opened with openbabel.\n Exiting...")
        exit()
     
    builder = solventBoxBuilder(solutexyz, solvent=solvent, slu_netcharge=slu_netcharge, cube_size=cube_size, charge_method=charge_method, 
                                slu_spinmult=slu_spinmult, outputFile=outputFile, srun_use=srun_use, 
                                gaussianexe=gaussianexe, gaussiandir=gaussiandir, amberhome=amberhome, 
                                closeness=closeness, solvent_off=solvent_off, solvent_frcmod=solvent_frcmod)
    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startboxgen(argumentList)