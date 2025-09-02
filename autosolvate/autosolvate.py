from openbabel import pybel
import getopt, sys, os
from openbabel import openbabel as ob
import subprocess 
import pkg_resources
from autosolvate.globs import keywords_avail, available_qm_programs, available_charge_methods
from autosolvate.resp_classes.resp_factory import resp_factory
from autosolvate.pubchem_api import PubChemAPI
from autosolvate.solute_info import Solute
from autosolvate.ionsFF import file_prep_for_ion


amber_solv_dict = {'water': [' ','TIP3PBOX '],
                   'methanol': ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
                   'chloroform': ['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
                   'nma': ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}

custom_solv_dict = {'acetonitrile':'ch3cn'}

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
            charge_method="resp", slu_spinmult=1, outputFile="", 
            srun_use=False, qm_program='gaussian', qm_exe=None, qm_dir=None,
            amberhome=None, closeness=0.8,
            solvent_off="", solvent_frcmod="",rundir=""):
        self.xyz = xyzfile
        self.solute = pybel.readfile('xyz', xyzfile).__next__()
        self.slu_netcharge = slu_netcharge
        self.slu_spinmult = slu_spinmult
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
        self.slu_pos = self.cube_size/2.0
        self.slv_count = 210*8 # 210
        self.pbcbox_size = self.cube_size+2
        self.outputFile=outputFile
        if self.outputFile == "":
            self.outputFile = self.solvent + "_solvated"
        self.srun_use=srun_use
        self.charge_method=charge_method
        self.qm_program = qm_program
        self.qm_exe=qm_exe
        self.qm_dir = qm_dir
        self.amberhome = amberhome
        self.rundir = rundir
        self.solvent_off=solvent_off
        self.solvent_frcmod=solvent_frcmod
        self.is_custom_solvent = False
        self.inputCheck()

    def inputCheck(self):
        if self.charge_method not in available_charge_methods:
            print("Error: the requested charge method: {}, is not supported yet".format(self.charge_method)) 
            exit(1)
        if self.slu_spinmult > 1:
            if self.charge_method != "resp":
                print("Error: solute spin multiplicity: ", self.slu_spinmult, " charge method: ", self.charge_method)
                print("Error: atomic charge fitting for open-shell system only works for resp charge method")
                print("Error: exiting...")
                exit(1)

        if self.charge_method == "resp":
            if self.qm_program not in available_qm_programs:
               print("Error! Selectd quantum chemistry package for RESP charge fittig is not supported:", qm_program)
               exit()
                
            elif self.qm_program == "gaussian":
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
                    if not os.path.exists(gspath):
                        print("Error! Gaussian executable path",gspath)
                        print("Does not exist! Exiting...")
                        exit()

            elif self.qm_program == "gamess":
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
        
        if os.path.exists(self.xyz):
            with open(self.xyz) as f:
                lines = f.readlines()
                atomnumber = int(lines[0].strip())
                if atomnumber == 1:
                    print('Warning this is an ion!')
                    self.ion = True
                else:
                    self.ion = False
        else:
            print("Error: solute xyz file does not exist")
            exit()       
                    

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
        obConversion.WriteFile(obmol,os.path.join(self.rundir,'solute.xyz.pdb'))
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
        print("Remeove CONNECT information from pdb file")
        self.removeConectFromPDB()

        if self.charge_method == "resp":
           myresp = resp_factory(pdbfile="solute.xyz.pdb", charge=self.slu_netcharge,
                                 spinmult=self.slu_spinmult, qm_program=self.qm_program,
                                 qm_exe=self.qm_exe, qm_dir=self.qm_dir,srun_use=self.srun_use,rundir=self.rundir)
           myresp.run()

        if self.charge_method == "bcc":
           print("AnteChamber: Generate mol2 with bcc charge.")
           cmd3="$AMBERHOME/bin/antechamber -i solute.xyz.pdb -fi pdb -o solute.mol2 -fo mol2 -c bcc -eq 2 -rn SLU"
           cmd3 += " -nc " + str(self.slu_netcharge) + " -m " + str(self.slu_spinmult)
           print(cmd3)
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
            print("tleap failed to generate solute.pdb") 
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
        self.getHeadTail()
        f = open("leap_add_solventbox.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n")
        f.write("loadoff " + self.solvent_off + "\n")
        f.write("loadamberparams " + self.solvent_frcmod + "\n")
        f.write("loadamberparams solute.frcmod\n")
        f.write("mol=loadmol2 solute.mol2\n")
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
        solvPrefix = custom_solv_dict[self.solvent]
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
        if self.solvent not in  custom_solv_dict.keys():
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
            packmol_inp.write("structure solute.pdb\n")
            packmol_inp.write("   number 1\n")
            packmol_inp.write("   fixed " + " " + str(self.slu_pos) + " "+ str(self.slu_pos) + " " + str(self.slu_pos) + " 0. 0. 0.\n")
            packmol_inp.write("   centerofmass\n")
            packmol_inp.write("   resnumbers 2 \n")
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
        solvPrefix = custom_solv_dict[self.solvent]
        solvent_frcmod = solvPrefix+'.frcmod'
        solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvent_frcmod))

        solvent_prep = solvPrefix+'.prep'
        solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvent_prep))
        f = open("leap_packmol_solvated.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
        f.write("loadamberparams " + solvent_frcmod_path + "\n")
        f.write("loadAmberPrep " + solvent_prep_path + "\n")
        f.write("loadamberparams solute.frcmod\n")
        f.write("loadoff solute.lib\n")
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
        elif self.solvent in custom_solv_dict.keys():
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
        if self.ion:
            ionFF = file_prep_for_ion(xyzfile=self.xyz, charge=self.slu_netcharge, solvent=self.solvent, outputFile=self.outputFile,
                                      cubic_size=self.cube_size,closeness=self.closeness, solvent_frcmod=self.solvent_frcmod,solvent_off=self.solvent_off)
            ionFF.build()
            print("The script has finished successfully")
        else:
            self.getSolutePDB()
            self.getFrcmod()
            self.createLib()
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
         -n, --solutename  initialize suggestion for box parameter for a given solute's name
         -m, --main  solute xyz file
         -s, --solvent  name of solvent (water, methanol, chloroform, nma)
         -o, --output  prefix of the output file names
         -c, --charge  formal charge of solute
         -u, --spinmultiplicity  spin multiplicity of solute
         -g, --chargemethod  name of charge fitting method (bcc, resp)
         -b, --cubesize  size of solvent cube in angstroms
         -r, --srunuse  option to run inside a slurm job
         -q, --qmprogram name of the quantum chemistry program (gaussian, gamess, terachem)
         -e, --qmexe  name of the quantum chemistry package executable used to generate electrostatic potential needed for RESP charge fitting
         -d, --qmdir  path to the quantum chemistry package
         -a, --amberhome  path to the AMBER molecular dynamics package root directory. Definition of the environment variable $AMBERHOME
         -t, --closeness  Solute-solvent closeness setting, for acetonitrile tolerance parameter in packmol in Å, for water, methanol, nma, chloroform the scaling factor in tleap, setting to 'automated' will automatically set this parameter based on solvent.
         -l, --solventoff  path to the custom solvent .off library file. Required if the user want to use some custom solvent other than the 5 solvents contained in AutoSolvate (TIP3P water, methanol, NMA, chloroform, MeCN)
         -p, --solventfrcmod  path to the custom solvent .frcmod file. Required if the user wants to use some custom solvent other than the 5 solvents contained in AutoSolvate. 
         -v, --validation  option to run validation step for given solute's name
         -h, --help  short usage description

    Returns
    -------
    None
        Generates the structure files and save as ```.pdb```. Generates the MD parameter-topology and coordinates files and saves as ```.prmtop``` and ```.inpcrd```
    """

    #print(argumentList)
    options = "hm:n:s:o:c:b:g:u:rq:e:d:a:t:l:p:vD:"
    long_options = ["help", "main","solutename", "solvent", "output", "charge", "cubesize", "chargemethod", "spinmultiplicity", "srunuse","qmprogram","qmexe", "qmdir", "amberhome", "closeness","solventoff","solventfrcmod", "validation", "runningdirectory"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    solutename = ""
    solutexyz=""
    solvent='water'
    slu_netcharge=0
    cube_size=54
    charge_method="bcc"
    slu_spinmult=1
    mult_given = 0
    mult_suggest = 0
    outputFile=""
    srun_use=False
    amberhome=None
    qmprogram="gaussian"
    qmexe=None
    qmdir=None
    closeness=0.8
    solvent_off=""
    solvent_frcmod=""
    rundir = ""
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
            print('  -q, --qmprogram            name of the quantum chemistry program')
            print('  -e, --qmexe                name of the quantum chemistry package executable')
            print('  -d, --qmdir                path to the quantum chemistry package')
            print('  -a, --amberhome            path to the AMBER molecular dynamics package root directory')
            print('  -t, --closeness            Solute-solvent closeness setting')
            print('  -l, --solventoff           path to the custom solvent .off library file')
            print('  -D, --rundir               running directory where temporary files are stored')
            print('  -p, --solventfrcmod        path to the custom solvent .frcmod file')
            print('  -n, --solutename           initialize suggested parameter for given solute')
            print('  -v, --validation           option to run validation step for input parameters')
            print('  -h, --help                 short usage description')
            exit()
        elif currentArgument in ("-m", "--main"):
            print ("Main/solutexyz", currentValue)
            solutexyz=str(currentValue) 
        elif currentArgument in ("-n", "--solutename"):
            if solutexyz == "":
                solutename=str(currentValue)
                sol=PubChemAPI(solutename)
                info=sol.get_info()
                solutexyz=str(info[3])
                slu_netcharge = info[2]
                solS=Solute(info[0], info[1], info[2], info[3])
                cube_size = solS.get_box_length()
                slu_spinmult = solS.get_spin_multiplicity()
                charge_method = solS.get_methods()[0]
            else:
                solS = Solute("", "", slu_netcharge, solutexyz)
                mol = next(pybel.readfile("xyz", solutexyz))
                total_electrons = sum(atom.atomicnum for atom in mol.atoms)
                slu_spinmult = (total_electrons - slu_netcharge) % 2 + 1
                if slu_spinmult > 1: charge_method = 'resp'
                cube_size = solS.get_box_length()    
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
        elif currentArgument in ("-q","--qmprogram"):
            print("Quantum chemistry program name:", currentValue)
            qmprogram = currentValue
        elif currentArgument in ("-e","--qmexe"):
            print("Quantum chemistry package executable name:", currentValue)
            qmexe = currentValue
        elif currentArgument in ("-d","--qmdir"):
            print("Quantum chemistry package directory:", currentValue)
            qmdir = currentValue
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
        elif currentArgument in ("-v", "--validation"):
            print('Validating...')
            solS = Solute("", "", slu_netcharge, solutexyz)
            mol = next(pybel.readfile("xyz", solutexyz))
            total_electrons = sum(atom.atomicnum for atom in mol.atoms)
            mult_suggest = ((total_electrons - slu_netcharge) % 2 + 1) % 2
            mult_given = slu_spinmult % 2
            if slu_spinmult == 0 or mult_suggest != mult_given:
                raise Exception("Incorrect solute spin multiplicity given, please double check your value or use suggestion function enabled by -n or --solutename")
            if slu_spinmult > 1 and charge_method != 'resp':
                raise Exception("Incorrect charge method given, please double check your value or use suggestion function enabled by -n or --solutename")
            if cube_size < solS.get_box_length():
                raise Exception("Solvent box is too small, please increase box length.")
            print('Passed')
        elif currentArgument in ('-D, --rundir '):
            print("Directory for Temporary files:",currentValue)
            rundir = currentValue

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
     
    builder = solventBoxBuilder(solutexyz, solvent=solvent, slu_netcharge=slu_netcharge,
                                cube_size=cube_size, charge_method=charge_method, 
                                slu_spinmult=slu_spinmult, outputFile=outputFile, srun_use=srun_use, 
                                qm_program=qmprogram, qm_exe=qmexe, qm_dir=qmdir,
                                amberhome=amberhome, closeness=closeness, 
                                solvent_off=solvent_off, solvent_frcmod=solvent_frcmod,rundir=rundir)
    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startboxgen(argumentList)
