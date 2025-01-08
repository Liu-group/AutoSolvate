import getopt, sys, os, subprocess


from .general_docker import GeneralDocker
from ..molecule import *
from ..utils import srun, compute_total_charge
from ..FFmetalcomplex import genFF


"""
class genFF():
    def __init__(self, 
                 filename:str, metal_charge:int, spinmult:int, chargefile:str, totalcharge:int,     # solute property
                 amberhome:str, cutoff:float, fakecharge:bool, mode:str,                            # AutoMCPB parameters
                 basisset:str, method:str, software:str, QMexe:str, maxcore:int, nprocs:int, opt:bool,  # QM calculation 
                 solvent:str, slv_count:int, solvent_off:str, solvent_frcmod:str,     # solvent property, ignored if using new boxgen
                 cubesize:float, closeness:float,                                     # solvation box generation, ignored if using new boxgen
                 folder:str, # the folder store the temporary output files
                 ):
"""

class AutoMCPBDocker(GeneralDocker):

    def __init__(self, 
                 basisset:str="DEF2-TZVP", method:str="B3LYP", software:str="orca", QMexe:str=None, maxcore:int=1024, nprocs:int=1, opt:bool=True,  # QM calculation
                 amberhome:str = "", cutoff:float =2.8, fakecharge:bool = False, mode:str = "A",
                 workfolder:str = WORKING_DIR, exeoutfile:str = "automcpb.log"
                 ) -> None:
        super(AutoMCPBDocker, self).__init__(
            executable = "python",
            workfolder = workfolder,
            exeoutfile = exeoutfile)
        
        self.logger.name = self.__class__.__name__

        # QM calculation
        self.basisset = basisset
        self.method = method
        self.software = software
        if not QMexe:
            QMexe = subprocess.check_output(["which", software]).decode().strip()
        self.QMexe = QMexe
        self.maxcore = maxcore
        self.nprocs = nprocs
        self.opt = opt

        # AutoMCPB parameters
        if not amberhome:
            amberhome = "$AMBERHOME/bin/"
        self.amberhome = amberhome
        self.cutoff = cutoff
        self.fakecharge = fakecharge
        self.mode = mode

        self.generated_prmtop = ""
        self.mcpb_folder = ""

    def check_system(self, mol: TransitionMetalComplex) -> bool:
        if not isinstance(mol, TransitionMetalComplex):
            self.logger.error("The input molecule is not a TransitionMetalComplex object.")
            return False
        if not mol.check_exist("pdb") and not mol.check_exist("xyz"):
            self.logger.error("The input molecule does not have a valid pdb or xyz file.")
            return False
        if mol.check_exist("prmtop"):
            self.logger.warning(f"The prmtop file of {mol.reference_name} already exists, will be overwritten.")
        return True
    
    def predict_output(self, mol: TransitionMetalComplex) -> str:
        self.mcpb_folder = mol.name + "_automcpb"
        self.mcpb_folder = os.path.join(self.workfolder, self.mcpb_folder)
        self.generated_prmtop = mol.name + "_dry.prmtop"
        self.generated_pdb = mol.name + "_dry.pdb"
        self.charge_info = mol.name + ".info"

    def generate_input(self, mol):
        pass 
    
    def generate_cmd(self, mol: TransitionMetalComplex) -> str:
        pass 

    def execute(self, mol: TransitionMetalComplex) -> bool:
        # override the parent class 'execute' method
        genff = genFF(
            # contained in TransitionMetalComplex
            filename = mol.name, 
            metal_charge = mol.metal_charge,
            spinmult = mol.spinmult,
            chargefile = mol.charge_file,
            totalcharge= mol.totalcharge,

            # AutoMCPB parameters
            amberhome=self.amberhome,
            cutoff=self.cutoff,
            fakecharge=self.fakecharge,
            mode=self.mode,

            # QM calculation
            basisset    = str(self.basisset)  ,
            method      = str(self.method  )  ,
            software    = str(self.software)  ,
            QMexe       = str(self.QMexe   )  ,
            maxcore     = str(self.maxcore )  ,
            nprocs      = str(self.nprocs  )  ,
            opt         = "Y" if (self.opt == True or self.opt == "Y") else "N",

            # solvent property, ignored
            solvent="water", slv_count=0, solvent_off="", solvent_frcmod="",
            cubesize=0, closeness=0,

            # folder
            folder = self.mcpb_folder
        )
        genff.build(use_boxgen_metal=False)

    def check_output(self, mol: TransitionMetalComplex) -> bool:
        if not os.path.exists(os.path.join(self.mcpb_folder, self.generated_prmtop)):
            self.logger.error(f"The prmtop file {self.generated_prmtop} is not generated.")
            raise FileNotFoundError(f"The prmtop file {self.generated_prmtop} is not generated.")
        
        if not os.path.exists(os.path.join(self.mcpb_folder, self.generated_pdb)):
            self.logger.error(f"The pdb file {self.generated_pdb} is not generated.")
            raise FileNotFoundError(f"The pdb file {self.generated_pdb} is not generated.")

        if not os.path.exists(os.path.join(self.mcpb_folder, self.charge_info)):
            self.logger.error(f"The legand charge file {self.charge_info} is not generated.")
            raise FileNotFoundError(f"The clegand charge file {self.charge_info} is not generated.")
        
        return True
    
    def process_output(self, mol: TransitionMetalComplex) -> None:
        self.logger.info(f"Set the result folder o {mol.name} as {self.mcpb_folder}.")
        setattr(mol, "origin_folder", self.mcpb_folder)
        self.logger.info(f"Find the prmtop file {self.generated_prmtop}.")
        setattr(mol, "prmtop", os.path.join(self.mcpb_folder, self.generated_prmtop))
        self.logger.info(f"Find the pdb file {self.generated_pdb}.")
        setattr(mol, "pdb", os.path.join(self.mcpb_folder, self.generated_pdb))
        self.logger.info(f"Find the charge file {self.charge_info}.")
        setattr(mol, "charge_file", os.path.join(self.mcpb_folder, self.charge_info))
        mol.update()

        if mol.totalcharge.upper() in ["DEFAULT", "AUTO"]:
            mol.totalcharge = compute_total_charge(mol.charge_file)
            self.logger.info(f"The total charge of {mol.name} is calculated as {mol.totalcharge}.")
        elif int(mol.totalcharge) != compute_total_charge(mol.charge_file):
            self.logger.warning(f"The total charge of {mol.name} is not consistent with the charge file.")
            self.logger.warning(f"Please check the charge file {mol.charge_file}.")
        mol.charge = mol.totalcharge
        mol.check_system()

    def run(self, mol: TransitionMetalComplex) -> bool:
        if not self.check_system(mol):
            return False
        self.predict_output(mol)
        self.generate_input(mol)
        self.generate_cmd(mol)
        self.execute(mol)
        if not self.check_output(mol):
            return False
        self.process_output(mol)
    

        
    