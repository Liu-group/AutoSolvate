import getopt, sys, os, subprocess
from typing import TextIO

from ..Common import * 
from ..molecule import *
from ..utils import srun
from .general_docker import GeneralDocker



class PackmolDocker(GeneralDocker):

    _PACKMOL_NLOOP = 200
    _PACKMOL_MAX_REMEDY_TRIALS = 10
    _PACKMOL_BOX_INCREMENT_ANG = 1.0

    def __init__(self, 
                 workfolder:            str = WORKING_DIR,
                 exeoutfile:            str = None,
        ) -> None:
        '''
        @NOTE: 
        1. I try to not let a function oriented class initialize any variables 
        2. Do not initialize a list as instance variable, because it will be shared by all instances 
           (python will only create one list for all instances)

        @TODO: 
        1. set up running directory
        '''
           
        super().__init__(
            executable = 'packmol',
            workfolder = workfolder,
            exeoutfile = sys.stdout)
        self.logger.name    = self.__class__.__name__

        self.packmolinp        = os.path.join(self.workfolder, "packmol.inp")
        self.packmolout        = os.path.join(self.workfolder, "packmol.out")
        self.outpdb         = None

    def check_system(self, mol: SolvatedSystem) -> None:
        if not isinstance(mol, SolvatedSystem):
            self.logger.critical("The input system is not a SolvatedSystem instance!")
            raise TypeError("The input system is not a SolvatedSystem instance!")
        self.logger.info("Input system: {}".format(mol.name))
        self.logger.info("tolerance: {}".format(mol.closeness))
        solutes = mol.solutes
        solvents = mol.solvents
        if len(solutes) == 0:
            self.logger.info("The input system does not contain any solutes!")
        elif len(solutes) > 1:
            self.logger.info("The input system contains more than one kind of solutes. Position of solutes will be randomly generated.")
        if len(solvents) == 0:
            self.logger.critical("The input system does not contain any solvents!")
            raise ValueError("The input system does not contain any solvents!")
        for slu in solutes:
            if isinstance(slu, (Molecule, TransitionMetalComplex)):
                self.logger.info("Solute: {}".format(slu.name))
                self.logger.info("\t PDB: {}".format(slu.pdb))
                self.logger.info("\t charge: {}".format(slu.charge))
                self.logger.info("\t multiplicity: {}".format(slu.multiplicity))
                self.logger.info("\t solute count: {}".format(slu.number))
            elif isinstance(slu, MoleculeComplex):
                self.logger.info("Solute: {}".format(slu.name))
                self.logger.info("\t PDB: {}".format(slu.pdb))
                self.logger.info("\t charge: {}".format(slu.netcharge))
                self.logger.info("\t multiplicity: {}".format(slu.multiplicity))
                self.logger.info("\t solute count: {}".format(slu.number))
            if slu.number > 1:
                self.logger.info("The solute {} has more than one molecule, the position of solvents will be randomly generated.".format(slu.name))
            elif slu.centered:
                self.logger.info("The solute {} is centered at the box.".format(slu.name))
            elif not slu.centered:
                self.logger.info("The solute {}'s position will be randomly generated.".format(slu.name))
        for slv in solvents:
            if isinstance(slv, SolventBox):
                self.logger.critical("The solvent should also be a molecule, not a pre-built solvent box!")
                raise TypeError("The solvent should also be a molecule, not a pre-built solvent box!")
            self.logger.info("Solvent: {}".format(slv.name))
            self.logger.info("\t PDB: {}".format(slv.pdb))
            self.logger.info("\t charge: {}".format(slv.charge))
            self.logger.info("\t multiplicity: {}".format(slv.multiplicity))
            self.logger.info("\t solvent count: {}".format(slv.number))

    def predict_output(self, mol: SolvatedSystem) -> None:
        inpname = os.path.join(self.workfolder, "{}_packmol.inp".format(mol.name))
        outname = os.path.join(self.workfolder, "{}_packmol.out".format(mol.name))
        outpdb  = os.path.join(self.workfolder, "{}.pdb".format(mol.name))
        for name, fmt in zip([inpname, outname, outpdb], ["packmol input", "packmol output", "pdb"]):
            self.logger.info("The {} file will be generated at {}".format(fmt, name))
            if os.path.exists(name):
                self.logger.info("Found a existing file with the same name: {}".format(name))
                self.logger.info("This file will be Overwritten!".format(name))
        self.packmolinp = inpname
        self.packmolout = outname
        self.outpdb     = outpdb
        self.output_files = [self.packmolinp, self.packmolout, self.outpdb]

    # generate input file
    def load_solute(self, doc: TextIO, box: SolvatedSystem) -> None: 
        r'''
        1. Already supports multiple solutes
        2. Already simplified

        @QUESTION
        1. does have to use 'solute.pdb' here? 
        2. check if {com} is correct
        @ANSWER
        1. No, just use the molecule's pdb file
        2. {com} is actually the rotation angle in degree. 
        '''
        solutes = box.solutes
        if len(solutes) == 0:
            doc.write('{}                                        \n'.format('# no solute'))
            return 
        solute = solutes[0]
        solute: Molecule
        center_pos = box.cubesize / 2.0 
        if len(solutes) == 1 and solute.number == 1:
            
            doc.write('{}                                        \n'.format('# add the solute'))
            # doc.write('{:<15} {}                                 \n'.format('structure', solute.pdb))
            if not os.path.exists(os.path.join(self.workfolder, os.path.basename(solute.pdb))):
                shutil.copy(solute.pdb, os.path.join(self.workfolder, os.path.basename(solute.pdb)))
                # this is actually a bug in packmol. The path to the pdb file cannot be too long. One have to copy the file to the working directory
            doc.write('{:<15} {}                                 \n'.format('structure', os.path.basename(solute.pdb)))
            doc.write('{:<15} {:<5}                              \n'.format('number',    1))
            doc.write('{:<10} {posx} {posy} {posz} {com} {com} {com}\n'.format('fixed', posx=center_pos[0], posy=center_pos[1], posz=center_pos[2], com=0.0))
            doc.write('{:<15}                                    \n'.format('centerofmass'))
            doc.write('{:<15} {:<5}                              \n'.format('resnumbers', '2'))
            doc.write('{:<15}                                    \n'.format('end structure'))
            doc.write('\n')
        else:
            for solute in solutes:
                doc.write('# add solute {}                           \n'.format(solute.name))
                # doc.write('{:<15} {}                                 \n'.format('structure', solute.pdb))
                if not os.path.exists(os.path.join(self.workfolder, os.path.basename(solute.pdb))):
                    shutil.copy(solute.pdb, os.path.join(self.workfolder, os.path.basename(solute.pdb)))
                doc.write('{:<15} {}                                 \n'.format('structure', os.path.basename(solute.pdb)))
                doc.write('{:<15} {:<5}                              \n'.format('number',solute.number))

                if solute.number == 1 and solute.centered:
                    doc.write('{:<10} {posx} {posy} {posz} {com} {com} {com}\n'.format('fixed', posx=center_pos[0], posy=center_pos[1], posz=center_pos[2], com=0.0))
                    doc.write('{:<15}                                    \n'.format('centerofmass'))
                else:
                    doc.write('{:<10} {xmin} {ymin} {zmin} {xmax} {ymax} {zmax}\n'.format('inside box', xmin='0.', ymin='0.', zmin='0.', xmax=box.cubesize[0], ymax=box.cubesize[1], zmax=box.cubesize[2]))
                doc.write('{:<15} {:<5}                              \n'.format('resnumbers', '2'))
                doc.write('{:<15}                                    \n'.format('end structure'))
                doc.write('\n')
        doc.write('\n')
        
    def load_solvent(self, doc: TextIO, box: SolvatedSystem) -> None: 
        solvents = box.solvents
        for solvent in solvents:
            solvent: Molecule
            doc.write('{} {}                                    \n'.format('# add solvent', solvent.name))
            # doc.write('{:<15} {}                                \n'.format('structure', solvent.pdb))
            if not os.path.exists(os.path.join(self.workfolder, os.path.basename(solvent.pdb))):
                shutil.copy(solvent.pdb, os.path.join(self.workfolder, os.path.basename(solvent.pdb)))
            doc.write('{:<15} {}                                \n'.format('structure', os.path.basename(solvent.pdb)))
            doc.write('{:<15} {:<5}                             \n'.format('number',    solvent.number))
            doc.write('{:<10} {xmin} {ymin} {zmin} {xmax} {ymax} {zmax}\n'.format('inside box', xmin='0.', ymin='0.', zmin='0.', xmax=box.cubesize[0], ymax=box.cubesize[1], zmax=box.cubesize[2]))
            doc.write('{:<15} {:<5}                             \n'.format('resnumbers', '2'))
            doc.write('{:<15}                                   \n'.format('end structure'))
            doc.write('\n')    

    def write_packmol_inp(self, box: SolvatedSystem) -> None: 
        r'''
        @TODO 
        1. finish this function 
        '''
        f = open(self.packmolinp, 'w') 
        f.write('{:<15} {:<3.2f}                \n'.format('tolerance', box.closeness))
        f.write('{:<15} {:<5}                   \n'.format('nloop', self._PACKMOL_NLOOP))
        f.write('{:<15} {:<5}                   \n'.format('filetype', 'pdb'))
        f.write('{:<15} {:<5}                   \n'.format('output', self.outpdb)) 
        self.load_solute(f, box)
        self.load_solvent(f, box)
        f.close()

    def _packmol_output_text(self) -> str:
        if not os.path.exists(self.packmolout):
            return ""
        with open(self.packmolout, 'r') as f:
            return f.read()

    def _packmol_started(self) -> bool:
        content = self._packmol_output_text()
        if not content:
            return False
        if "PACKMOL - Packing optimization" in content:
            return True
        if "Starting GENCAN loop" in content:
            return True
        return False

    def _expand_box(self, mol: SolvatedSystem, delta_ang: float) -> None:
        mol.cubesize = mol.cubesize + delta_ang
        self.logger.warning(
            "Packmol remedy: expanded box by %.2f angstrom. New box size: %.2f %.2f %.2f",
            delta_ang,
            mol.cubesize[0],
            mol.cubesize[1],
            mol.cubesize[2],
        )

    def generate_input(self, mol: SolvatedSystem) -> None:
        self.write_packmol_inp(mol)
    @srun()
    def generate_cmd(self) -> str:
        cmd = "{} < {} > {}".format(self.executable, self.packmolinp, self.packmolout)
        return cmd
    
    # execute of packmol is a little different with other softwares:
    # - packmol prefers running in workfolder
    # - input/output paths should be relative and co-located
    def execute(self, cmd: str) -> None:
        tmp_packmolinp = os.path.basename(self.packmolinp)
        tmp_packmolout = os.path.basename(self.packmolout)
        cmd = "{} < {} > {}".format(self.executable, tmp_packmolinp, tmp_packmolout)

        # Reuse GeneralDocker's milestone prints + stderr routing, while
        # keeping packmol's relative-path constraint.
        prev = getattr(self, "exeoutfile", None)
        try:
            self.exeoutfile = None
            super().execute(cmd)
        finally:
            self.exeoutfile = prev

    
    # check if the output file is generated        
    def check_output(self):
        if not os.path.exists(self.packmolout):
            self.logger.critical("The packmol output file {} is not found!".format(self.packmolout))
            raise FileNotFoundError("The packmol output file {} is not found!".format(self.packmolout))
        with open(self.packmolout, 'r') as f:
            content = f.read()
        if content.find("                   Success!") == -1:
            # this is for falling back to the generated pdb even if packmol did not fully converge, which is a common issue. The generated pdb may not be perfect, but it is better than nothing.
            self.logger.warning("Packmol failed to generate the system pdb!")
            # self.logger.critical("Please check the packmol output file for details.")
            raise RuntimeError("Packmol failed to generate the system pdb!")
        if not os.path.exists(self.outpdb):
            self.logger.critical("The system pdb {} is not found!".format(self.outpdb))
            raise FileNotFoundError("The system pdb {} is not found!".format(self.outpdb))
        self.logger.info("Packmol successfully generated the system pdb: {}".format(self.outpdb))
    
    def process_output(self, mol:SolvatedSystem):
        self.logger.info("Adding TER section in the packmol generated pdb {} ...".format(self.outpdb))
        data   = open(self.outpdb,'r').readlines()
        output = open(self.outpdb,'w') 
        output : TextIO
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
        self.logger.info("The system pdb has been successfully processed.".format(self.outpdb))
        mol.pdb = self.outpdb
        mol.update()

    def run(self, mol: SolvatedSystem) -> None:
        self.logger.name = self.__class__.__name__
        self.check_system(mol)
        self.predict_output(mol)

        packmol_ok = False
        for remedy_trial in range(self._PACKMOL_MAX_REMEDY_TRIALS + 1):
            self.generate_input(mol)
            cmd = self.generate_cmd()
            self.execute(cmd)

            try:
                self.check_output()
                packmol_ok = True
                break
            except RuntimeError as err:
                if not self._packmol_started():
                    self.logger.critical(
                        "Packmol did not start correctly. This is not a packmol convergence failure; no remedy retry is applied."
                    )
                    self.logger.critical(f"Error string: {str(err)}")
                    raise err

                if remedy_trial < self._PACKMOL_MAX_REMEDY_TRIALS:
                    self._expand_box(mol, self._PACKMOL_BOX_INCREMENT_ANG)
                    self.logger.warning(
                        "Packmol attempt %d/%d failed after start; retrying with expanded box.",
                        remedy_trial + 1,
                        self._PACKMOL_MAX_REMEDY_TRIALS,
                    )
                    continue

                if os.path.exists(self.outpdb):
                    self.logger.warning(
                        "Packmol failed strict success criteria after %d remedy attempts, but generated pdb exists: %s. "
                        "Proceeding with this suboptimal structure.",
                        self._PACKMOL_MAX_REMEDY_TRIALS,
                        self.outpdb,
                    )
                    packmol_ok = True
                    break
                raise err

        if not packmol_ok:
            raise RuntimeError("Packmol failed and no usable pdb was generated.")

        self.process_output(mol)

