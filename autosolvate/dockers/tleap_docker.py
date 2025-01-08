
import getopt, sys, os, subprocess
from typing import TextIO

from ..Common import * 
from ..molecule import *
from ..utils import *
from .general_docker import GeneralDocker


class TleapDocker(GeneralDocker):

    _ff14SB         = 'leaprc.protein.ff14SB'       # Source leaprc file for ff14SB protein force field
    _gaff           = 'leaprc.gaff'                 # Source leaprc file for gaff force field
    _water_tip3p    = 'leaprc.water.tip3p'          # Source leaprc file for tip3p water model
    _ionslm_126_opc = "frcmod.ionslm_126_opc"        #Source frcmod file for ionslm_126_opc force field

    def __init__(self, 
                 workfolder:            str = WORKING_DIR,
                 exeoutfile:            str = None,
                 boxsize:               float = None,
    ) -> None:
        super().__init__(
            executable = 'tleap',
            workfolder = workfolder,
            exeoutfile = sys.stdout)
        self.logger.name    = self.__class__.__name__
        self.boxsize        = boxsize if isinstance(boxsize, Iterable) else [boxsize, boxsize, boxsize]
        self.boxsize        = None if not self.boxsize[0] else self.boxsize

        self.leapinp        = os.path.join(self.workfolder, "leap.inp")
        self.leapout        = os.path.join(self.workfolder, "leap.out")

        self.outpdb         = None
        self.outlib         = None
        self.outmol2        = None
        self.outprmtop      = None
        self.outinpcrd      = None

    ##### check_system method for different systems
    def check_system_molecule(self, mol:Molecule):
        self.logger.info("Checking system {:s}...".format(mol.name))
        if mol.amber_solvent:
            self.logger.info("This is a predefined AMBER solvent.")
            if mol.name == "water":
                self.logger.info(f"Using TIP3P water model with predefined residue name {mol.residue_name}")
            else:
                self.logger.info("Solvent prep file: {:s}".format(mol.prep))
                self.logger.info("Solvent frcmod file: {:s}".format(mol.frcmod))
            return 
            
        suffixs = ["mol2", "lib", "prep", "off"]
        suffixf = [0, 0, 0, 0]
        for i, suffix in enumerate(suffixs):
            if mol.check_exist(suffix):
                self.logger.info("\t System {}: {:s}".format(suffix, getattr(mol, suffix)))
                suffixf[i] = 1
            else:
                self.logger.debug("\t System {} not found".format(suffix, getattr(mol, suffix)))
        if sum(suffixf) > 1:
            self.logger.info("More than one structure files were found, parameters may be overwritten.")
        suffix = "frcmod"
        if mol.check_exist("frcmod"):
            self.logger.info("\t System {}: {:s}".format(suffix, getattr(mol, suffix)))
        else:
            self.logger.warning("System {} not found".format(suffix, getattr(mol, suffix)))
            self.logger.warning("Tleap may failed to generate AMBER parameter files!")
            self.logger.warning("Please check the tleap output file {}".format(self.leapout))
        if sum(suffixf) <= 0:
            self.logger.critical("None of the {} files were found!".format(suffixs))
            raise RuntimeError("None of the {} files were found!".format(suffixs))
        if self.boxsize is not None:
            self.logger.info("Add a box with shape {}".format(self.boxsize))
        
    def check_system_multimolecule(self, mol:MoleculeComplex):
        self.logger.info("Checking system {:s}...".format(mol.name))
        self.logger.info("This is a molecule complex.")
        self.logger.info("Number of molecules in this system: {}".format(len(mol.fragmols)))
        self.logger.info("New molecules in this system: {}".format([m.name for m in mol.newmolecules]))
        for m in mol.newmolecules:
            m:Molecule
            self.logger.info("New molecule: {}"     .format(m.name))
            self.logger.info("\t PDB: {}"           .format(m.pdb))
            self.logger.info("\t charge: {}"        .format(m.charge))
            self.logger.info("\t multiplicity: {}"  .format(m.multiplicity))
            self.logger.info("\t residue name: {}"  .format(m.residue_name))
            self.check_system_molecule(m)
        
    def check_system_solvatedsystem(self, mol:SolvatedSystem):
        self.logger.info("Checking system {:s}...".format(mol.name))
        if isinstance(mol.cubesize, float):
            self.logger.info("Cube size: {:.2f}".format(mol.cubesize))
            self.boxsize = [mol.cubesize, mol.cubesize, mol.cubesize]
        elif isinstance(mol.cubesize, Iterable):
            self.logger.info("Box size: {:.2f} {:.2f} {:.2f}".format(*mol.cubesize))
            self.boxsize = mol.cubesize
        self.logger.info("Solutes in this system: {}".format([slu.name for slu in mol.solutes]))
        for slu in mol.solutes:
            slu:System
            self.logger.info("Checking solute {:s}...".format(slu.name))
            self.logger.info("Number of {}: {}".format(slu.name, slu.number))
            if isinstance(slu, Molecule):
                self.check_system_molecule(slu)
            elif isinstance(slu, MoleculeComplex):
                self.check_system_multimolecule(slu)
        for slv in mol.solvents:
            self.logger.info("Checking solvent {:s}...".format(slv.name))
            if isinstance(slv, SolventBox):
                if os.path.exists(slv.solventbox):
                    self.logger.info("Pre-built solvent box: {}".format(slv.solventbox))
                elif slv.amber_solvent:
                    self.logger.info("Use amber solvent box: {}".format(slv.solventbox))
                else:
                    self.logger.critical("Solvent box not found: {}".format(slv.solventbox))
                    raise RuntimeError("Solvent box not found: {}".format(slv.solventbox))
                if os.path.exists(slv.frcmod):
                    self.logger.info("Pre-built solvent frcmod: {}".format(slv.frcmod))
                elif slv.amber_solvent:
                    self.logger.info("Use amber solvent frcmod: {}".format(slv.frcmod))
                else:
                    self.logger.critical("Solvent frcmod not found: {}".format(slv.frcmod))
                    raise RuntimeError("Solvent frcmod not found: {}".format(slv.frcmod))
            elif isinstance(slv, Molecule): 
                self.logger.info("Custom solvent detected")
                self.logger.info("Number of {}: {}".format(slv.name, slv.number))
                self.check_system_molecule(slv)
                if not mol.check_exist("pdb"):
                    self.logger.critical("PDB file not found for this solvated system!")
                    self.logger.critical("Please use Packmol to generate a PDB file for this system before running tleap!")
                    raise RuntimeError("PDB file not found for this solvated system!")
        
    def check_system(self, mol:System):
        if isinstance(mol, Molecule):
            self.check_system_molecule(mol)
        elif isinstance(mol, MoleculeComplex):
            mol.tleap_pre_process()
            self.check_system_multimolecule(mol)
        elif isinstance(mol, SolvatedSystem):
            self.check_system_solvatedsystem(mol)
        elif isinstance(mol, TransitionMetalComplex):
            pass             # the prmtop and inpcrd files are already generated by AutoMCPB
        else:
            self.logger.critical("Unknown system type!" + str(type(mol)))
            raise TypeError("Unknown system type!" + str(type(mol)))

    ##### predict_output
    def get_output_name(self, mol:System, fmt:str):
        return mol.reference_name + "." + fmt

    def predict_output(self, mol:System):
        for fmt in ["pdb", "lib", "mol2", "prmtop", "inpcrd"]:
            outname = self.get_output_name(mol, fmt)
            self.logger.debug("The {} file will be generated at {}".format(fmt, outname))
            if os.path.exists(outname):
                self.logger.info("Found a existing file with the same name: {}".format(outname))
                self.logger.info("This file will be Overwritten!".format(outname))
            setattr(self, "out"+fmt, outname)
        self.leapinp = os.path.join(self.workfolder, "leap_{}.cmd".format(mol.name))
        self.logger.debug("The tleap input file will be generated at {}".format(self.leapinp))
        self.leapout = os.path.join(self.workfolder, "leap_{}.log".format(mol.name))
        self.logger.debug("The tleap output file will be generated at {}".format(self.leapout))
    
    ##### generate_input method for different systems
    def load_forcefield(self, doc: TextIO) -> None:
        doc.write('{:<20}  {:<20}       \n'.format('source', self._ff14SB)) 
        doc.write('{:<20}  {:<20}       \n'.format('source', self._gaff))
        doc.write('{:<20}  {:<20}       \n'.format('source', self._water_tip3p))

    def load_mol(self,  
                 doc:           TextIO, 
                 mol:           System, 
                 check:         bool = False
    ) -> None: 
        r'''
        @QUESTION:
        1. what is check mol for? 
        '''
        if isinstance(mol, Molecule):
            if mol.check_exist("frcmod"): 
                doc.write('{:<20}  {}       \n'.format('loadamberparams', mol.frcmod))       
            if mol.check_exist("lib"): 
                doc.write('{:<20}  {}       \n'.format('loadoff', mol.lib)) 
            if mol.check_exist("off"):
                doc.write('{:<20}  {}       \n'.format('loadoff', mol.off))
            if mol.check_exist("prep"): 
                doc.write('{:<20}  {}       \n'.format('loadamberprep', mol.prep))
            if mol.check_exist("mol2"):  
                doc.write('{:<5} {:<15} {}  \n'.format(mol.residue_name, '= loadmol2', mol.mol2))
            if check: 
                doc.write('{:<20}  {}       \n'.format('check', mol.residue_name))
        elif isinstance(mol, MoleculeComplex):
            for m in mol.newmolecules:
                m:Molecule
                if m.check_exist("frcmod"): 
                    doc.write('{:<20}  {}       \n'.format('loadamberparams', m.frcmod))
                if m.check_exist("lib"): 
                    doc.write('{:<20}  {}       \n'.format('loadoff', m.lib))
            doc.write('{:<5} {:<15} {}          \n'.format(mol.residue_name, '= loadpdb', mol.pdb))
        elif isinstance(mol, TransitionMetalComplex):
            mol.write_tleap_load_command(doc)
            if mol.check_exist("pdb"):
                doc.write('{:<5} {:<15} {}          \n'.format(mol.residue_name, '= loadpdb', mol.pdb))
                mol.write_tleap_add_bond_command(doc, mol.pdb, mol.residue_name)

  
    def load_solvent_box(self,
                        doc:           TextIO,
                        mol:           SolventBox,
                        check:         bool = False
    ) -> None:
        if mol.amber_solvent and not mol.name == "water":
            doc.write('{:<5} {}  \n'.format("loadamberparams", mol.frcmod))
            return 
        
        if mol.check_exist("off"):
            doc.write('{:<5} {}  \n'.format("loadoff", mol.off))
        if mol.check_exist("lib"):
            doc.write('{:<5} {}  \n'.format("loadoff", mol.lib))
        if mol.check_exist("frcmod"):
            doc.write('{:<5} {}  \n'.format("loadamberparams", mol.frcmod))
        if check:
            doc.write('{:<5} {}  \n'.format("check", mol.name))

    def load_head_tail(self, doc:TextIO, mol:Molecule) -> None:
        head, tail = getHeadTail(mol.mol2)
        if head == 0 and tail == 0:
            return 
        doc.write('{} {} {} {}.{}.{} \n'.format("set", mol.residue_name, "head", mol.residue_name, 1, head))
        doc.write('{} {} {} {}.{}.{} \n'.format("set", mol.residue_name, "tail", mol.residue_name, 1, tail))


    def write_tleap_in_transition_metal_complex(self, mol:TransitionMetalComplex) -> None:
        self.logger.info("Tleap input file: {}".format(self.leapinp))
        f = open(self.leapinp, 'w')
        self.load_forcefield(f)
        self.load_mol(f, mol)
        f.write('{} {} {}        \n'.format('savepdb', mol.residue_name, self.outpdb))
        f.write('{} {} {} {}     \n'.format('saveamberparm', mol.residue_name, self.outprmtop, self.outinpcrd)) 
        f.write('{}              \n'.format('quit'))
        f.close()

    def write_tleap_in_molecule(self, mol:Molecule) -> None:
        r'''
        @QUESTION 
        1. Is check mol necessary?, what is it for? 
        2. Does self count as an object? in @dispatch(object) 

        @TODO 
        1. check if 20 pt is enough for all strings 
        '''
        self.logger.info("Tleap input file: {}".format(self.leapinp))
        f = open(self.leapinp, 'w')
        self.load_forcefield(f)
        self.load_mol(f, mol) 
        if mol.check_exist("mol2"):
            self.load_head_tail(f, mol) 
        if self.boxsize is not None:
            f.write('set %s box {%.2f,%.2f,%.2f} \n' % (mol.residue_name, self.boxsize[0], self.boxsize[1], self.boxsize[2]))
        f.write('{} {} {} {}     \n'.format('savemol2', mol.residue_name, self.outmol2, 1))
        f.write('{} {} {}        \n'.format('saveoff', mol.residue_name, self.outlib))
        f.write('{} {} {}        \n'.format('savepdb', mol.residue_name, self.outpdb))
        f.write('{} {} {} {}     \n'.format('saveamberparm', mol.residue_name, self.outprmtop, self.outinpcrd)) 
        f.write('{}              \n'.format('quit'))
        f.close() 
    
    def write_tleap_in_add_solvent(self, mol:SolvatedSystem) -> None:
        self.logger.info("Tleap input file: {}".format(self.leapinp))
        f = open(self.leapinp, 'w')
        self.load_forcefield(f)
        slu: Molecule   = mol.solutes[0]
        slv: SolventBox = mol.solvents[0]
        slu_pos = mol.cubesize / 2
        self.load_mol(f, slu)
        self.load_solvent_box(f, slv)
        mol_tleap_name = slu.residue_name
        f.write('{} {} {} {} {} {}\n'.format(
            'solvatebox',
            slu.residue_name, 
            slv.box_name, 
            slu_pos[0],
            "iso",
            mol.closeness,))
        if mol.netcharge != 0:
            ion = "Cl-" if mol.netcharge > 0 else "Na+"
            f.write('{} {} {} {} \n'.format('addions',slu.residue_name,ion,0))
        f.write('{} {} {}        \n'.format('savepdb', mol_tleap_name, self.outpdb))
        f.write('{} {} {} {}     \n'.format('saveamberparm', mol_tleap_name, self.outprmtop, self.outinpcrd)) 
        f.write('{}              \n'.format('quit'))
        f.close()

    def write_tleap_in_custom_solvated(self, mol:SolvatedSystem) -> None:
        self.logger.info("Tleap input file: {}".format(self.leapinp))
        f = open(self.leapinp, 'w')
        self.load_forcefield(f)
        for slu in mol.solutes:
            self.load_mol(f, slu)
        for slv in mol.solvents:
            self.load_mol(f, slv)
        mol_tleap_name = "SYS"
        f.write('{} = {} {}      \n'.format(mol_tleap_name, "loadpdb", mol.pdb))
        for slu in mol.solutes:
            if isinstance(slu, TransitionMetalComplex):
                slu.write_tleap_add_bond_command(f, mol.pdb, mol_tleap_name)
        if mol.netcharge != 0:
            ion = "Cl-" if mol.netcharge > 0 else "Na+"
            f.write('{} {} {} {} \n'.format('addions', mol_tleap_name, ion, 0))
        pbcbox_size = mol.cubesize + 2
        f.write("set %s box {%.2f,%.2f,%.2f}\n" % (mol_tleap_name, pbcbox_size[0], pbcbox_size[1], pbcbox_size[2]))
        f.write('check {}        \n'.format(mol_tleap_name))
        f.write('{} {} {}        \n'.format('savepdb', mol_tleap_name, self.outpdb))
        f.write('{} {} {} {}     \n'.format('saveamberparm', mol_tleap_name, self.outprmtop, self.outinpcrd)) 
        f.write('{}              \n'.format('quit'))
        f.close()

    def generate_input(self, mol:System):
        if isinstance(mol, Molecule) or isinstance(mol, MoleculeComplex):
            self.write_tleap_in_molecule(mol)
        elif isinstance(mol, TransitionMetalComplex):
            self.write_tleap_in_transition_metal_complex(mol)
        elif isinstance(mol, SolvatedSystem):
            if isinstance(mol.solvent, SolventBox):
                self.write_tleap_in_add_solvent(mol)
            else:
                self.write_tleap_in_custom_solvated(mol)

    ##### generate_cmd method
    @srun()
    def generate_cmd(self, mol:System) -> str:
        cmd = 'tleap -s -f {} > {}'.format(self.leapinp, self.leapout)
        return cmd
    
    ##### check_output
    def check_tleap_output(self):
        with open(self.leapout, "r") as f:
            content = f.read()
        if content.find("Exiting LEaP: Errors = 0") == -1:
            eline = ""
            for line in content.splitlines():
                if line.find("FATAL") != -1:
                    eline = line
            raise AssertionError("Error in running tleap. " + eline)

    def check_output(self, mol:System):
        success = True
        for fmt, outfile in zip(["prmtop", "inpcrd"], [self.outprmtop, self.outinpcrd]):
            if os.path.exists(outfile):
                self.logger.info("Successfully generated target {} file: {}".format(fmt, outfile))
                success = True
            else:
                self.logger.critical("Failed to generate target {} file: {}".format(fmt, outfile))
                success = False
        if not success:
            raise RuntimeError("Failed to generate target files!")
        self.check_tleap_output()

    ##### process_output
    def process_output(self, mol:System):
        mol.prmtop = self.outprmtop
        mol.inpcrd = self.outinpcrd
        mol.lib    = self.outlib
        if self.outpdb:
            mol.pdb = self.outpdb
        mol.update()

    def run(self, mol:System):
        self.logger.name = self.__class__.__name__
        self.check_system(mol)
        self.predict_output(mol)
        self.generate_input(mol)
        cmd = self.generate_cmd(mol)
        self.execute(cmd)
        self.check_output(mol)
        self.process_output(mol)
