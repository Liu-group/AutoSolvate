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
# Update 2024-07-12:
#   1. Add a MixtureBuilder class to generate mixed solventbox for both single solute and molecular pairs
#   2. Now accept pre-generated prep, lib, and off files if the frcmod file is provided. 
#   3. Mixed solvent with amber defined solvent such as TIP3P water, methanol are enabled. 
#--------------------------------------------------------------------------------------------------#
import getopt, sys, os
import subprocess
from typing import List, Tuple, Iterable
import json
import logging
import inspect
import argparse

from .molecule import *
from .dockers import *
from .utils import *

from autosolvate.autosolvate import *

class MulticomponentParamsBuilder():
    """
    Create amber parameter files for a single xyz or pdb file with multiple separate fragments.

    Warning: If you want to create the forcefield for transition metal complexes please use the ```boxgen_metal``` module instead of ```boxgen_multicomponent```.

    Parameters
    ----------
    xyzfile : str
        structure file name, can be any structural files that openbabel recognizes.
    name : array_like, Optional. default: the base name of the provided structure file.
        
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
    """
    Build a solvent box for a single molecule complex as the solute and single solvent

    Parameters
    ----------
    xyzfile : str
        structure file of the molecular complex, can be any type within ["xyz", "pdb", "mol2"]. "prep", "lib", "off" are not supported for molecular complex.
    slu_charge : int or dict, Optional, default: 0
        Charge of the solute. if complex has charged fragments, provide a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as neutral.
    slu_spinmult : int or dict, Optional, default: 1
        Spin multiplicity of the solute. if complex has non-singlet fragments, provide a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as singlet.
    charge_method : str, Optional, default: "resp"
        name of charge fitting method (bcc, resp)
    slu_count : int, Optional, default: 1
        number of the solute in the system. Not recommanded to set this parameter. Be cautious about this as the solute may have more than 1 fragments.
    solvent : str, Optional, default: "water"
        name of the solvent. Predefined solvents include ["water", "methanol", "chloroform", "nma"].
    solvent_frcmod : str, Optional, default: ""
        path to the frcmod file of the solvent. Required when user have the solvent forcefield parameters prepared.
    solvent_off : str, Optional, default: ""
        path to the off file of the solvent. Required when user have the solvent forcefield parameters prepared.
    solvent_box_name : str, Optional, default: "SLVBOX"
        name of the solvent box
    slv_generate : bool, Optional, default: False
        whether to generate the solvent forcefield parameters with GAFF. If True, the solvent forcefield will be generated with GAFF, where 'slv_xyz' parameter will be needed. If False, the solvent will be treated as a predefined solvent in AMBER or user should provide the frcmod & prep files.
    slv_xyz : str, Optional, default: ""
        path to the xyz file of the solvent. Required when user want to generate the solvent forcefield parameters with GAFF.
    slv_count : int, Optional, default: 210*8
        number of the solvent in the system.
    cube_size : int, Optional, default: 54
        size of solvent cube in angstroms
    closeness : float, Optional, default: 0.8
        Solute-solvent closeness setting, corresponding to the tolerance parameter in packmol in Å,
    folder : str, Optional, default: current working directory
        the directory where the files are generated
    outputFile : str, Optional, default: ""
        prefix of the output .pdb, .inpcrd and .prmtop files
    kwargs : dict
        Other arguments that need to be included in the solute. Remained for future development.    
    """
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
        """
        Build the solvated system with molecule complex as the solute and single solvent 

        Parameters
        ----------
        None

        Returns
        -------
        None
        
        """
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
    prefix : str, Optional, default: None
        prefix of the output file names. Default will be <solute_name>_<solvent_name_1>_...-<solvent_name_n>
    """
    def __init__(self, folder = WORKING_DIR, cube_size = 54, closeness = 2.0, charge_method = "bcc", prefix = None):
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
        self.systemprefix = prefix

    def add_transition_metal_complex_solute(self, xyzfile:str, number = 1, 
            metal_charge = 0, spinmult = 1, total_charge = "default", chargefile = "", 
            qm_kwargs:dict = {}, mcpb_kwargs:dict = {}, **kwargs):
        """
        add a transition metal complex as the solute. e.g. a metalloenzyme

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any type within ["xyz", "pdb", "mol2"]. "prep", "lib", "off" are not supported for transition metal complex.
        number : int, Optional, default: 1
            number of the solute in the system. Not recommanded to set this parameter. Be cautious about this as the solute may have more than 1 fragments.
        metal_charge : int, Optional, default: 0
            charge of the metal atom, must be provided.
        spinmult : int, Optional, default: 1
            spin multiplicity of the compound, must be provided.
        total_charge : int, Optional, default: "default"
            total charge of the complex. If not provided, the total charge will be calculated from the charge of the metal atom and the chargefile.
        chargefile : str, Optional, default: ""
            path to the charge file of the legands. if not provided, all legand will be considered as neutral.
        qm_kwargs : dict, Optional, default: {}
            additional arguments for the quantum chemistry calculation of the metal complex. Default parameters:
            {
                "method": "b3lyp",
                "basisset": "DEF2-TZVP",
                "software": "orca",
                "QMexe": <where the orca program installed>,
                "maxcore": 1024,
                "nprocs": 1,
                "opt": True,
                }
        mcpb_kwargs : dict, Optional, default: {}
            additional arguments for the mcpb.py calculation of the metal complex. Default parameters:
            {
                "amberhome": $AMBERHOME/bin/,
                "cutoff": 2.8,
                "fakecharge": False,
                "mode": "A",
            }
        """
        
        tmc = TransitionMetalComplex(xyzfile, folder = self.folder, metal_charge = metal_charge, multiplicity = spinmult, totalcharge = total_charge, legand_charge_file = chargefile)
        dockerparams = qm_kwargs.copy()
        dockerparams.update(mcpb_kwargs)
        dockerparams["workfolder"] = self.folder
        self.logger.info(dockerparams)
        mcpb_docker = AutoMCPBDocker(**dockerparams)
        mcpb_docker.run(tmc)
        tmc.number = number
        self.solutes.append(tmc)

    def add_complex_solute(self, xyzfile:str, fragment_charge = 0, fragment_spinmult = 1, number = 1, **kwargs):
        """
        add a molecular complex as the solute. e.g. a electron transfer donor-acceptor pair

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any type within ["xyz", "pdb", "mol2"]. "prep", "lib", "off" are not supported for molecular complex.
        fragment_charge : dict | array_like, Optional, default: 0
            Charge for each fragment. A dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as neutral.
        fragment_spinmult : dict | array_like, Optional, default: 1
            Multiplicity for each fragment. A dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as singlet.
        number : int, Optional, default: 1
            number of the solute in the system. Not recommanded to set this parameter. Be cautious about this as the solute may have more than 1 fragments.
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
            additional files needed for the solute, including "mol2", "frcmod", "lib", "prep", and "off". Will support "itp", "top" in the future.
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
              ("off" in kwargs and os.path.isfile(kwargs["off"])) or \
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
            structure file name, can be any type within ["xyz", "pdb", "mol2", "prep", "off"]. Can be ignored if the solvent is predefined in AMBER.
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
            if "off" in kwargs and os.path.isfile(kwargs["off"]):
                self_solvent.off = kwargs["off"]
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
        """
        Start to build the mixed solvent box. No parameters are needed.

        Parameters
        ----------
        None

        Returns
        -------
        None

        """
        if not self.systemprefix:
            system_name = "-".join([m.name for m in self.solutes] + [m.name for m in self.solvents])
        else:
            system_name = self.systemprefix
        self.solutes :List[Molecule] 
        self.solvents:List[Molecule]
        solute_numbers  = [m.number for m in self.solutes ]
        solvent_numbers = [m.number for m in self.solvents]
        system = SolvatedSystem(system_name, solute = self.solutes, solvent = self.solvents,
                                cubesize=self.boxsize, closeness=self.closeness, 
                                solute_number = solute_numbers, solvent_number = solvent_numbers,
                                folder = self.folder)
        for docker in self.custom_solvation:
            docker.run(system)




def startmulticomponent_fromdata(data:dict):
    """
    Start the multicomponent solvation process from a python dictionary.

    Parameters
    ----------
    data : dict
        dictionary containing the input parameters. Usually generated from a json file.
    """

    parser = InputParser()
    parser.read_dict(data)
    parser.parse()
    data = parser.data

    data["folder"] = WORKING_DIR
    json.dump(data, open(os.path.join(data["folder"], "autosolvate_input_full.json"), "w"), indent=4)

    signature = inspect.signature(MixtureBuilder.__init__)
    function_params = signature.parameters
    filtered_data = {k: v for k, v in data.items() if k in function_params}
    builder = MixtureBuilder(**filtered_data)
    for solute in data["solutes"]:
        if solute["__TYPE__"] == "molecule":
            builder.add_solute(**solute)
        elif solute["__TYPE__"] == "complex":
            builder.add_complex_solute(**solute)
        elif solute["__TYPE__"] == "transition_metal_complex":
            if "qm_kwargs" in data:
                solute["qm_kwargs"] = data["qm_kwargs"]
            if "mcpb_kwargs" in data:
                solute["mcpb_kwargs"] = data["mcpb_kwargs"]
            builder.add_transition_metal_complex_solute(**solute)
    for solvent in data["solvents"]:
        builder.add_solvent(**solvent)
    builder.build()

def startmulticomponent_fromfile(file:str):
    """
    Start the multicomponent solvation process from a json file.

    Parameters
    ----------
    file : str
        json file containing the input parameters. 
    """
    with open(file, "r") as f:
        data = json.load(f)
    startmulticomponent_fromdata(data)

def create_parser_multicomponent():
    parser = argparse.ArgumentParser(
        description='Add solvent box to a given solute and generate related force field parameters.',
        epilog="suggest usage: autosolvate multicomponent -f <JSON path> \nif an input file is provided, all command line options will be ignored. \nIf using command line as the traditional way, it will only generate a single solute with single solvent. \nThis is a legacy feature, designed solely for the compatibility with the older version. It is not recommended for further use."
    )

    parser.add_argument('-f', '--file',            type=str,  help='json file containing the input parameters. Will ignore all other options if provided. Required when using multiple solvents')
    parser.add_argument('-m', '--main',            type=str,  default='',      help='solute xyz file')
    parser.add_argument('-o', '--output',          type=str,  default='',      help='prefix of the output file names')
    parser.add_argument('-c', '--charge',          type=int,  default=0,       help='formal charge of solute')
    parser.add_argument('-u', '--spinmultiplicity',type=int,  default=1,       help='spin multiplicity of solute')
    parser.add_argument('-s', '--solvent',         type=str,  default='water', help='solvent xyz files, Will use single solvent if provided.')
    parser.add_argument('-g', '--chargemethod',    type=str,  default='bcc',   help='name of charge fitting method (bcc, resp)')
    parser.add_argument('-b', '--cubesize',        type=float,default=54.0,    help='size of solvent cube in angstroms')
    parser.add_argument('-t', '--closeness',       type=float,default=2.0,     help='solute-solvent closeness setting. Automation is not possible for mixed solvent')    
    parser.add_argument('-r', '--srunuse',         action='store_true',        help='option to run inside a slurm job')
    parser.add_argument('-e', '--gaussianexe',     type=str,                   help='name of the Gaussian quantum chemistry package executable')
    parser.add_argument('-d', '--gaussiandir',     type=str,                   help='path to the Gaussian package')
    parser.add_argument('-a', '--amberhome',       type=str,                   help='path to the AMBER molecular dynamics package root directory')
    return parser

def startmulticomponent(args):
    r"""
    Wrap function that parses command line options for autosolvate multicomponent module,
    generate solvent box and related force field parameters.
    suggested usage: autosolvate multicomponent -f <JSON path>
    
    Command Line Options
    --------------------
        -f, --file 
            json file containing the input parameters, Required when using multiple solvents. Will ignore all other options if provided.
        -m, --main  
            solute xyz file, An Legacy feature, designed for the compatibility with the older version. It is not recommended for further use.
        -o, --output  
            prefix of the output file names
        -c, --charge  
            formal charge of solute
        -u, --spinmultiplicity  
            spin multiplicity of solute
        -s, --solvent  
            solvent xyz files, Will use single solvent if provided. Not available for using multiple solvents
        -g, --chargemethod  
            name of charge fitting method (bcc, resp)
        -b, --cubesize  
            size of solvent cube in angstroms
        -t, --closeness  
            Solute-solvent closeness setting. Default 2.0 Å for mixed solvent. For acetonitrile tolerance parameter in packmol in Å, for water, methanol, nma, chloroform the scaling factor in tleap, setting to 'automated' will automatically set this parameter based on solvent.
        -r, --srunuse  
            option to run inside a slurm job
        -e, --gaussianexe  
            name of the Gaussian quantum chemistry package executable used to generate electrostatic potential needed for RESP charge fitting
        -d, --gaussiandir  
            path to the Gaussian package
        -a, --amberhome  
            path to the AMBER molecular dynamics package root directory. Definition of the environment variable $AMBERHOME
        -h, --help  
            short usage description

    Returns
    -------
        Generates the structure files and save as ```.pdb```. Generates the MD parameter-topology and coordinates files and saves as ```.prmtop``` and ```.inpcrd```
    """
    #print(argumentList)
    parser = create_parser_multicomponent()
    args = parser.parse_args(args)
    cmd_dict = vars(args)
    if "file" in cmd_dict and cmd_dict["file"] and os.path.exists(cmd_dict["file"]):
        if not os.path.exists(cmd_dict["file"]):
            raise FileNotFoundError(f"File {cmd_dict['file']} not found")
        with open(cmd_dict["file"], "r") as f:
            data = json.load(f)
    else:
        cmd_dict.pop("file")
        data = convert_cmd_to_dict(cmd_dict)
    global WORKING_DIR
    WORKING_DIR = os.getcwd()
    startmulticomponent_fromdata(data)
    

if __name__ == "__main__":
    startmulticomponent(sys.argv[1:])
    
