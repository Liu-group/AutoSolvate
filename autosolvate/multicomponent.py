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
import json
import logging
import inspect

from .molecule import *
from .dockers import *
from .utils import *

from autosolvate.autosolvate import *

class MulticomponentParamsBuilder():
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
    def __init__(self, 
                 xyzfile: str, name="", residue_name="SLU", charge=0, spinmult=1, 
                 charge_method="resp", folder = WORKING_DIR, **kwargs): 

        
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
        if solvent in AMBER_SOLVENTBOX_DICT and not slv_generate: # predefined amber solvent
            self_solvent = AMBER_SOLVENTBOX_DICT[solvent]
            self_solvent.folder = self.folder
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
    """
    Create amber parameter files for a single solute with mixed solvents

    Use 'add_solute' to add solute and 'add_solvent' to add solvent.

    Parameters
    ----------
    folder : str, Optional, default: current working directory
        working directory
    cube_size : int, Optional, default: 54
        size of solvent cube in angstroms
    closeness : float, Optional, default: 2.0
        Solute-solvent closeness setting, corresponding to the tolerance parameter in packmol in Å, 
    charge_method : str, Optional, default: "bcc"
        name of charge fitting method (bcc, resp)
    """
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
        self.logger = logging.getLogger(name = self.__class__.__name__)
        self.output_handler             = logging.FileHandler(filename = "autosolvate.log", mode = "a", encoding="utf-8")
        self.output_formater            = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
        self.output_handler.setFormatter(self.output_formater)
        if len(self.logger.handlers) == 0:
            self.logger.addHandler(self.output_handler)


    def add_complex_solute(self, xyzfile:str, fragment_charge = 0, fragment_spinmult = 1, number = 1, **kwargs):
        """
        add a molecular complex as the solute 

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any type within ["xyz", "pdb", "mol2"]. frcmod are not supported for molecular complex.
        fragment_charge : dict | array_like, Optional, default: 0
            Charge for each fragment. A dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as neutral.
        fragment_spinmult : dict | array_like, Optional, default: 1
            Multiplicity for each fragment. A dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as singlet.
        number : int, Optional, default: 1
            number of the solute in the system. Be cautious about this as the solute may have more than 1 fragments.
        **kwargs : dict
            Other arguments that need to be included in the solute. Remained for future development.
        """
        molecule = MoleculeComplex(xyzfile, charges=fragment_charge, multiplicities=fragment_spinmult, folder = self.folder)
        for fragment in molecule.newmolecules:
            for docker in self.single_molecule_pipeline:
                docker.run(fragment)
        TleapDocker(workfolder = self.folder).run(molecule)
        molecule.number = number
        self.solutes.append(molecule)
        
    def add_solute(self, xyzfile:str, name="", residue_name="SLU", charge=0, spinmult=1, number = 1, **kwargs):
        """
        add a solute molecule

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any type within ["xyz", "pdb", "mol2"]
        name : str, Optional, default: the base name of the provided structure file
            name of the solute
        residue_name : str, Optional, default: "SLU"
            residue name of the solute. Note if an mol2 or prep file is provided, the residue name will be read from the file.
        charge : int, Optional, default: 0
            charge of the solute
        spinmult : int, Optional, default: 1
            spin multiplicity of the solute
        number : int, Optional, default: 1
            number of the solute in the system.
        **kwargs : dict
            additional files needed for the solute, including "mol2", "frcmod", "lib", "prep".
            If the user want to skip the antechamber and leap steps, the user need to provide the "mol2" and "frcmod" files by adding the following arguments:
            mol2 : str  
                the path of the mol2 file of the solute
            frcmod : str
                the path of the frcmod file of the solute
        """

        if ("mol2" in kwargs and os.path.isfile(kwargs["mol2"])) and \
           ("frcmod" in kwargs and os.path.isfile(kwargs["frcmod"])):
            molecule = Molecule(xyzfile, charge=charge, multiplicity=spinmult, folder = self.folder, name = name, residue_name=residue_name)
            molecule.mol2 = kwargs["mol2"]
            molecule.frcmod = kwargs["frcmod"]
            molecule.get_residue_name()
            molecule.update()
        else:
            molecule = Molecule(xyzfile, charge=charge, multiplicity=spinmult, folder = self.folder, name = name, residue_name=residue_name)
            for docker in self.single_molecule_pipeline:
                docker.run(molecule)
        molecule.number = number
        self.solutes.append(molecule)

    def get_solvent_type(self, xyzfile = "", name = "", **kwargs):
        solvent_type = "generate"
        if name in AMBER_SOLVENT_DICT:
            if (("prep" in kwargs and os.path.isfile(kwargs["prep"])) or \
                ("mol2" in kwargs and os.path.isfile(kwargs["mol2"]))) and \
                ("frcmod" in kwargs and os.path.isfile(kwargs["frcmod"])):
                solvent_type = "custom"
            else:
                solvent_type = "amber"
        elif name in custom_solv_dict:
            if (("prep" in kwargs and os.path.isfile(kwargs["prep"])) or \
                ("mol2" in kwargs and os.path.isfile(kwargs["mol2"]))) and \
                ("frcmod" in kwargs and os.path.isfile(kwargs["frcmod"])):
                solvent_type = "custom"
            else:
                solvent_type = "autosolvate_custom"
        elif (("prep" in kwargs and os.path.isfile(kwargs["prep"])) or \
              ("mol2" in kwargs and os.path.isfile(kwargs["mol2"]))) and \
              ("frcmod" in kwargs and os.path.isfile(kwargs["frcmod"])):
            solvent_type = "custom"
        elif os.path.exists(xyzfile):
            solvent_type = "generate"
        else:
            raise ValueError("Solvent not found")
        return solvent_type

    def add_solvent(self, xyzfile:str = "", name="", residue_name="SLV", charge=0, spinmult=1, number = 210*8, **kwargs):

        """
        add a type of solvent

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any type within ["xyz", "pdb", "mol2", "prep"]. Can be ignored if the solvent is predefined in AMBER.
        name : str
            name of the solvent. Predefined solvents include ["water", "methanol", "chloroform", "nma"].
        residue_name : str, Optional, default: "SLV"
            residue name of the solvent. Will be ignored if the solvent is predefined in AMBER. This argument will be ignored if a mol2 or prep file is provided together with a frcmod. 
        charge : int, Optional, default: 0
            charge of the solvent
        spinmult : int, Optional, default: 1
            spin multiplicity of the solvent
        number : int, Optional, default: 210*8
            number of the solvent in the system.
        **kwargs : dict
            additional files needed for the solvent, including "mol2", "frcmod", "lib", "prep", will support "itp", "top" in the future.
            If the user want to skip the antechamber and leap steps, the user need to provide the ["mol2" or "prep"] and "frcmod" files by adding the following arguments:
            mol2 : str  
                the path of the mol2 file of the solvent
            frcmod : str
                the path of the frcmod file of the solvent
        """

        solvent_type = self.get_solvent_type(xyzfile, name, **kwargs)
        if solvent_type == "amber":
            self.logger.info(f"Adding predefined solvent {name}")
            self_solvent = AMBER_SOLVENT_DICT[name]
            self_solvent.folder = self.folder
            self_solvent.generate_pdb()
        elif solvent_type == "autosolvate_custom":
            self.logger.info(f"Adding autosolvate custom solvent {name}")
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
        elif solvent_type == "custom":
            self.logger.info(f"Adding user provided custom solvent {name}")
            self_solvent = Molecule(xyzfile, charge=charge, multiplicity=spinmult, folder = self.folder, name = name, residue_name=residue_name)
            if "mol2" in kwargs and os.path.isfile(kwargs["mol2"]):
                self_solvent.mol2 = kwargs["mol2"]
            if "frcmod" in kwargs and os.path.isfile(kwargs["frcmod"]):
                self_solvent.frcmod = kwargs["frcmod"]
            if "lib" in kwargs and os.path.isfile(kwargs["lib"]):
                self_solvent.lib = kwargs["lib"]
            if "prep" in kwargs and os.path.isfile(kwargs["prep"]):
                self_solvent.prep = kwargs["prep"]
            self_solvent.get_residue_name()
            self_solvent.update()
        elif solvent_type == "generate":
            self.logger.info(f"Adding solvent {name} whose forcefield parameters will be generated with GAFF")
            self_solvent = Molecule(xyzfile, charge=charge, multiplicity=spinmult, folder = self.folder, name = name, residue_name=residue_name)
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

def startmulticomponent_fromfile(jsonpath:str):
    with open(jsonpath, "r") as f:
        data = json.load(f)
    data["folder"] = os.getcwd()
    signature = inspect.signature(MixtureBuilder.__init__)
    function_params = signature.parameters
    filtered_data = {k: v for k, v in data.items() if k in function_params}
    builder = MixtureBuilder(**filtered_data)

    if "solute" in data:
        data["solute"] = add_missing_xyzfile_keyword(data["solute"])
        if check_multicomponent(data["solute"]["xyzfile"]):
            builder.add_complex_solute(**data["solute"])
        else:
            builder.add_solute(**data["solute"])
    if "solvent" in data:
        data["solvent"] = add_missing_xyzfile_keyword(data["solvent"])
        builder.add_solvent(**data["solvent"])
    if "solutes" in data and type(data["solutes"]) == list:
        for i in range(len(data["solutes"])):
            data["solutes"][i] = add_missing_xyzfile_keyword(data["solutes"][i])
            if check_multicomponent(data["solutes"][i]["xyzfile"]):
                builder.add_complex_solute(**data["solutes"][i])
            else:
                builder.add_solute(**data["solutes"][i])
    if "solvents" in data and type(data["solvents"]) == list:
        for i in range(len(data["solvents"])):
            data["solvents"][i] = add_missing_xyzfile_keyword(data["solvents"][i])
            builder.add_solvent(**data["solvents"][i])
    builder.build()

    

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
    startmulticomponent(sys.argv[1:])
