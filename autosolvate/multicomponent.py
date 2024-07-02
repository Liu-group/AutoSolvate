#--------------------------------------------------------------------------------------------------#
# multicomponent.py. 
# Description: 
#   This module can handle structure files containing multiple molecules.
#   1. Generates the lib and frcmod file for each separate molecule when it is not h2o or amino acids.
#   2. Generates the lib, prmtop and inpcrd files of this whole structure.
# Update 2022-02-04:
#   1. Can generate lib, mol2 and prmtop file for an xyz with mutiple fragments with charges
#   2. Can create solvent box for xyz file with multiple fragments.
#   3. Able to rearrange shuffled pdb files into ordered form that antechamber can process.
# author: Fangning Ren (2022-02-04) 
# path: autosolvate/multicomponent.py
#--------------------------------------------------------------------------------------------------#
import getopt, sys, os
import subprocess
from typing import List, Tuple, Iterable

from .molecule import *
from .dockers import *
from .utils import *

from autosolvate.autosolvate import *

class MulticomponentParamsBuilder():
    def __init__(self, 
                 xyzfile: str, name="", residue_name="SLU", charge=0, spinmult=1, 
                 charge_method="resp", folder = WORKING_DIR, **kwargs): 
        """
        Create amber parameter files for a single xyz or pdb file with multiple separate fragments

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any structural files that openbabel recognizes.
        name : array_like, Optional. default: the base name of the provided structure file.
            Not used
        residue_name : array_like, Optional. default: Residue name provided in pdb file or assigned as UAA, UAB, UAC, etc.
            Residue names for each fragments. A list of strings of three capital letters. Its length should equal to the number of fragments in xyzfile. If this parameter is not given, the residues will be assigned by "U" plus "AB","AC",..."AZ", "BA"...
        charge : dict | array_like, Optional. default: 0
            Charge for each fragment. A list of integer with its length equal to the number of fragments, or a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as neutral. 
        spinmult : dict | array_like, Optional. default: 0
            Multiplicity for each fragment. A list of integer with its length equal to the number of fragments, or a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as singlet. 
        outputFile : str, Optional, default='water_solvated'
            Filename-prefix for outputfiles
        pre_optimize_fragments : bool, Optional, default: False
            do geometry optimization with MMFF94 forcefield in OpebBabel before running antechamber
        srun_use : bool, Optional, default='False
            Run all commands with a srun prefix
        gaussianexe : str, Optional, default: g16
            name of the Gaussian executeble
        gaussiandir : str, Optional, default: $GAUSSIANDIR
            path of Gaussian
        amberhome : str, Optional, default: $AMBERHOME
            path of amber
        deletefiles : bool, Optional, default: False
            Delete all temporary files except the .prmtop and .inpcrd file of the pdb file provided.
        """
        
        self.folder = folder
        self.mol = MoleculeComplex(xyzfile, name=name, residue_name=residue_name, charges=charge, multiplicities=spinmult, folder = self.folder)
        self.charge_method = charge_method

        self.single_molecule_pipeline = [
            AntechamberDocker(charge_method = self.charge_method, workfolder = self.folder),
            ParmchkDocker(workfolder=self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.complex_pipeline = [TleapDocker(workfolder=self.folder)]

    def build(self):
        for m in self.mol.newmolecules:
            if self.charge_method == "resp":
                build_resp_terachem(m, folder = self.folder)
            else:
                for docker in self.single_molecule_pipeline:
                    docker.run(m)
        for docker in self.complex_pipeline:
            docker.run(self.mol)

       
class MulticomponentSolventBoxBuilder():
    def __init__(self, 
                 xyzfile:str, slu_charge=0, slu_spinmult=1, charge_method="resp", slu_count = 1,
                 solvent = "water", solvent_frcmod = "", solvent_off = "", solvent_box_name = "SLVBOX",
                 slv_generate = False, slv_xyz = "", slv_count = 210*8,
                 cube_size = 54, closeness = 0.8, folder = WORKING_DIR, outputFile = "",
                 **kwargs):

        self.kwargs = kwargs

        if not outputFile:
            outputFile = solvent + "_solvated"
        self.folder = folder
        if "slu_netcharge" in kwargs and isinstance(kwargs["slu_netcharge"], dict) and slu_charge == 0:
            slu_charge = kwargs["slu_netcharge"]
        self.solute = MoleculeComplex(xyzfile, slu_charge, slu_spinmult, folder = self.folder)
        self.charge_method = charge_method

        self.solvent = self.get_solvent(solvent, slv_xyz, solvent_frcmod, solvent_off, slv_generate, slv_count, solvent_box_name)
        self.system = SolvatedSystem(solvent + "_solvated", solute = self.solute, solvent=self.solvent,
                                     cubesize=cube_size, closeness=closeness, solute_number=slu_count, solvent_number=slv_count,
                                     folder = self.folder)

        self.single_molecule_pipeline = [
            AntechamberDocker(charge_method = self.charge_method, workfolder = self.folder),
            ParmchkDocker(workfolder=self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.complex_pipeline = [TleapDocker(workfolder=self.folder)]
        self.solvent_pipeline = [
            AntechamberDocker("bcc", "mol2", workfolder = self.folder),
            ParmchkDocker("frcmod", workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.prebuilt_solvation = [
            TleapDocker(workfolder = self.folder)
        ]
        self.custom_solvation = [
            PackmolDocker(workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]
    
    def get_solvent(self, solvent:str, slv_xyz:str = "", solvent_frcmod:str = "", solvent_off:str = "", slv_generate:bool = False, slv_count:int = 210*8, solvent_box_name:str = "SLVBOX"):
        if solvent in AMBER_SOLVENT_DICT and not slv_generate:
            self_solvent = AMBER_SOLVENT_DICT[solvent]
        elif solvent in custom_solv_dict:
            # solvent data prepared by autosolvate
            solvPrefix = custom_solv_dict[solvent]
            solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".frcmod"))
            solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".prep"))
            solvent_pdb_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".pdb"))
            self_solvent = Molecule(solvent_pdb_path, 0, 1, solvent, residue_name = custom_solv_residue_name[solvent], folder = self.folder)
            self_solvent.frcmod = solvent_frcmod_path
            self_solvent.prep   = solvent_prep_path
            self_solvent.number = slv_count
        elif os.path.exists(slv_xyz) and slv_generate:
            # generate solvent box
            self_solvent = Molecule(slv_xyz, folder = self.folder)
            self_solvent.number = slv_count
        else:
            raise ValueError("Solvent not found")
        return self_solvent

    def build(self):
        for m in self.solute.newmolecules:
            if self.charge_method == "resp":
                build_resp_terachem(m, folder = self.folder)
            else:
                for docker in self.single_molecule_pipeline:
                    docker.run(m)
        for docker in self.complex_pipeline:
            docker.run(self.solute)

        if isinstance(self.solvent, Molecule):
            for docker in self.solvent_pipeline:
                docker.run(self.solvent)
        if isinstance(self.solvent, SolventBox):
            self.prebuilt_solvation[0].run(self.system)
        elif isinstance(self.solvent, Molecule):
            self.custom_solvation[0].run(self.system)
            self.custom_solvation[1].run(self.system)


class MixtureBuilder():
    def __init__(self, folder = WORKING_DIR, cube_size = 54, closeness = 2.0, charge_method = "bcc"):
        self.solutes = []
        self.solvents = []
        self.folder = folder
        self.boxsize = [cube_size, cube_size, cube_size]
        self.closeness = closeness
        self.charge_method = charge_method
        self.single_molecule_pipeline = [
            AntechamberDocker(charge_method = self.charge_method, workfolder = self.folder),
            ParmchkDocker(workfolder=self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.custom_solvation = [
            PackmolDocker(workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]

    def add_solute(self, xyzfile:str, name="", residue_name="", charge=0, spinmult=1, number = 1, **kwargs):

        # use the first three letters of the xyzfile name as the residue name if not provided 
        if not residue_name: 
            residue_name = try_ones_best_to_get_residue_name(xyzfile, name)
        molecule = Molecule(xyzfile, charge=charge, multiplicity=spinmult, folder = self.folder, name = name, residue_name=residue_name)

        if "mol2" in kwargs and os.path.isfile(kwargs["mol2"]):
            molecule.mol2 = kwargs["mol2"]
        if "frcmod" in kwargs and os.path.isfile(kwargs["frcmod"]):
            molecule.frcmod = kwargs["frcmod"]
        if "lib" in kwargs and os.path.isfile(kwargs["lib"]):
            molecule.lib = kwargs["lib"]
        molecule.update()

        if not molecule.check_exist("frcmod") or not (molecule.check_exist("mol2") or molecule.check_exist("lib")):
            for docker in self.single_molecule_pipeline:
                docker.run(molecule)
        molecule.number = number
        self.solutes.append(molecule)

    def add_solvent(self, xyzfile:str = "", name="", residue_name="", charge=0, spinmult=1, number = 210*8, slv_generate = False, **kwargs):
        if name in AMBER_SOLVENT_DICT and not slv_generate: # predefined amber solvent
            self_solvent = AMBER_SOLVENT_DICT[name]
            self_solvent.folder = self.folder
            self_solvent.generate_pdb()
        elif name in custom_solv_dict:                      # custom solvent in autosolvate
            # solvent data prepared by autosolvate
            solvPrefix = custom_solv_dict[name]
            solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".frcmod"))
            solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".prep"))
            solvent_pdb_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvPrefix+".pdb"))
            self_solvent = Molecule(solvent_pdb_path, 0, 1, name, residue_name = custom_solv_residue_name[name], folder = self.folder)
            self_solvent.frcmod = solvent_frcmod_path
            self_solvent.prep   = solvent_prep_path
        elif os.path.exists(xyzfile) and slv_generate:      # user defined solvent
            # generate solvent box
            residue_name = try_ones_best_to_get_residue_name(xyzfile, name)
            self_solvent = Molecule(xyzfile, charge=charge, multiplicity=spinmult, folder = self.folder, name = name, residue_name=residue_name)
            if "mol2" in kwargs and os.path.isfile(kwargs["mol2"]):
                self_solvent.mol2 = kwargs["mol2"]
            if "frcmod" in kwargs and os.path.isfile(kwargs["frcmod"]):
                self_solvent.frcmod = kwargs["frcmod"]
            if "lib" in kwargs and os.path.isfile(kwargs["lib"]):
                self_solvent.lib = kwargs["lib"]
            self_solvent.update()

            if not self_solvent.check_exist("frcmod") or not (self_solvent.check_exist("mol2") or self_solvent.check_exist("lib")):
                for docker in self.single_molecule_pipeline:
                    docker.run(self_solvent)
        else:
            raise ValueError("Solvent not found")
        self_solvent.number = number
        for solvent in self.solvents:
            if solvent.name == self_solvent.name:
                raise ValueError(f"Solvent {solvent.name} already exists")
        self.solvents.append(self_solvent)

    def build(self):
        system_name = "-".join([m.name for m in self.solutes] + [m.name for m in self.solvents])
        solute_numbers = [m.number for m in self.solutes]
        solvent_numbers = [m.number for m in self.solvents]
        '''
        @DEBUG 
        comment left by Patrick Jun 12 2024.  
        the old way to define 'system_name' will cause bug, please have the person who wrote this to fix it.
        '''
        '''
        @DEBUG
        bug fixed by Fangning Ren on July 1 2024
        '''
        system = SolvatedSystem(system_name, solute = self.solutes, solvent = self.solvents,
                                cubesize=self.boxsize, closeness=self.closeness, 
                                solute_number = solute_numbers, solvent_number = solvent_numbers,
                                folder = self.folder)
        for docker in self.custom_solvation:
            docker.run(system)
            
def startmulticomponent(argumentList):
    r"""
    Wrap function that parses command line options for autosolvate multicomponent module,
    generate solvent box and related force field parameters.
    
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
            print('Usage: autosolvate boxgen_multicomponent [OPTIONS]')
            print('  -m, --main                 solute xyz file')
            print('  -s, --solvent_lists        solvent xyz files, separated by slash. (ex. water.pdb/acetonitrile.pdb)')
            # print('  -o, --output               prefix of the output file names')                                     # not currently supporting this option
            print('  -c, --charge               formal charge of solute')            
            print('  -u, --spinmultiplicity     spin multiplicity of solute')
            print('  -g, --chargemethod         name of charge fitting method (bcc, resp)')                             # currently only support bcc 
            print('  -b, --cubesize             size of solvent cube in angstroms') 
            # print('  -r, --srunuse              option to run inside a slurm job')                                    # not currently supporting this option         
            # print('  -e, --gaussianexe          name of the Gaussian quantum chemistry package executable')           # not currently supporting gaussian
            # print('  -d, --gaussiandir          path to the Gaussian package')                                        # not currently supporting gaussian
            # print('  -a, --amberhome            path to the AMBER molecular dynamics package root directory')         # not currently supporting amberhome 
            print('  -t, --closeness            Solute-solvent closeness setting')                  
            # print('  -l, --solventoff           path to the custom solvent .off library file')                        # not currently supporting this option 
            # print('  -p, --solventfrcmod        path to the custom solvent .frcmod file')                             # not currently supporting this option     
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

    global WORKING_DIR
    WORKING_DIR = os.getcwd()
    
    builder = MixtureBuilder()
    builder.add_solute(solutexyz, charge=slu_netcharge, spinmult=slu_spinmult) 
    for s in solvent.split("/"):
        builder.add_solvent(s)
    builder.build()


if __name__ == "__main__":
    inst = MulticomponentParamsBuilder("PAHs.pdb", deletefiles=True)
    inst.build()
