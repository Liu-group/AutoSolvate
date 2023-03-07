from .molecule import *
from .dockers import *
from .utils import *


class AmberParamsBuilder(object):
    """ 
    This class handles the Amber parameter creation for one single molecule.
    1. Generate standard pdb
    2. AnteChamber or Gaussian charge fitting
    3. Tleap create Lib
    """
    def __init__(self, xyzfile:str, name = "", resname = "", charge = 0, spinmult = 1, 
                 charge_method="resp", folder = WORKING_DIR, **kwargs):
            
            self.folder = folder
            self.mol = Molecule(xyzfile, charge, spinmult, name = name, residue_name=resname, folder = self.folder)
            self.charge_method = charge_method

            self.bcc_pipeline = [
                AntechamberDocker("bcc", "mol2", workfolder = self.folder),
                ParmchkDocker("frcmod", workfolder = self.folder),
                TleapDocker(workfolder = self.folder)
            ]

            self.resp_pipeline = [] # 还没整

    def build(self, **kwargs):
        for docker in self.bcc_pipeline:
            docker:GeneralDocker
            docker.run(self.mol, **kwargs)

class solventBoxBuilder(object):
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
    def __init__(self, 
                 xyzfile:str, slu_netcharge=0, slu_spinmult=1, charge_method="resp", slu_count = 1,
                 solvent = "water", solvent_frcmod = "", solvent_off = "", solvent_box_name = "SLVBOX",
                 slv_generate = False, slv_xyz = "", slv_count = 210*8,
                 cube_size = 54, closeness = 0.8, 
                 **kwargs):
        
        self.kwargs = kwargs

        self.folder = os.path.splitext(xyzfile)[0] + "_solvated"

        self.solute = Molecule(xyzfile, slu_netcharge, slu_spinmult, folder = self.folder)
        self.solute.number = slu_count

        if solvent in AMBER_SOLVENT_NAMES:
            for slv in AMBER_SOLVENT_LIST:
                if slv.name == solvent:
                    self.solvent = slv
                    break
        elif os.path.exists(solvent_frcmod) and os.path.exists(solvent_off):
            self.solvent = SolventBox(solvent, solvent_off, solvent_frcmod, folder = self.folder, box_name=solvent_box_name)
        elif os.path.exists(slv_xyz) and slv_generate == True:
            self.solvent = Molecule(slv_xyz, folder = self.folder)
            self.solvent.number = slv_count
        else:
            raise ValueError("Solvent not found")
                # autosolvate solvent?

        self.system = SolvatedSystem(os.path.splitext(xyzfile)[0] + "_solvated", solute = self.solute, solvent=self.solvent,
                                     cubesize=cube_size, closeness=closeness, solute_number=slu_count, solvent_number=slv_count,
                                     folder = self.folder)
        self.solute_bcc_pipeline = [
            AntechamberDocker("bcc", "mol2", workfolder = self.folder),
            ParmchkDocker("frcmod", workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.solvent_pipeline = [
            AntechamberDocker("bcc", "mol2", workfolder = self.folder),
            ParmchkDocker("frcmod", workfolder = self.folder),
            TleapDocker(workfolder = self.folder)
        ]
        self.solvation_pipeline = [
            TleapDocker(workfolder = self.folder)
        ]
    
    def build(self):
        for docker in self.solute_bcc_pipeline:
            docker:GeneralDocker
            docker.run(self.solute)
        if isinstance(self.solvent, Molecule):
            for docker in self.solvent_pipeline:
                docker.run(self.solvent)
        self.solvation_pipeline[0].run(self.system)
    