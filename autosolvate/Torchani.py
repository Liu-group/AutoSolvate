from tempfile import TemporaryDirectory
from openbabel import pybel
from ssl import ALERT_DESCRIPTION_DECOMPRESSION_FAILURE
import sys, os
from openbabel import openbabel as ob
import subprocess
from types import CellType
import torchani
import ase
from ase.io import read, write
from ase.optimize import FIRE, BFGS
from ase import units
import subprocess
import pkg_resources
from ase import units
from ase.md.nptberendsen import NPTBerendsen
from ase.md.nvtberendsen import NVTBerendsen
from ase.md.langevin import Langevin
from ase.calculators.checkpoint import Checkpoint
from ase.io.trajectory import Trajectory


Number_A = 6.02214076 * 10**23

elements_dic = {'H':1,'HE':2,'LI':3,'BE': 4,'B':5,'C':6, 'N':7, 'O':8, 'F':9,'NE':10,
                 'NA': 10,'MG':12,'AL':13,'SI':14,'P':15,'S':16, 'CL':17, 'AR':18,'K':19,'CA':20}

weight_dic ={'water': 18.01528 } ##g/mol

density_dic = {'water' : 997 } ##kg/m3

def calculateSolvnumber(solvPrefix,volume):  ###Volume should be m3
    denisty = density_dic[solvPrefix]
    weight = weight_dic[solvPrefix]
    mass = volume * denisty  ## kg
    Mol = mass*1000/weight
    NumberOfSolv = Mol * Number_A
    return int(NumberOfSolv) 

class MLFF(object):
    def __init__(self, xyz, solvent='water.pdb', slu_netcharge=0, cube_size=20,closeness=0.8, tempi=0,temp0=90, dt=2, nstlim_heating=10000):
        self.xyz = xyz
        self.solvent = solvent
        self.closeness = closeness
        self.cube_size = cube_size
        self.slu_pos = self.cube_size/2.0
        self.pbcbox_size = self.cube_size+2
        self.solvPrefix = self.solvent.split('.pdb')[0]
        self.volumn = (self.cube_size * 10**(-10))**3 ### m
        self.slv_count = calculateSolvnumber(self.solvPrefix, self.volumn)
        self.tempi = tempi
        self.temp0 = temp0
        self.dt = dt
        self.nstlim_heating = nstlim_heating


    def packSLUSLV(self):

        solute_xyzfile = read(self.xyz)
       # pdbfilename = solute_xyzfile.split('.xyz')[0]+'.pdb'
        write('solute.pdb',solute_xyzfile)

        solvent_pdb = self.solvent
        output_pdb = self.solvPrefix + "_solvated.packmol.pdb"
    
    #   solvent_pdb_origin = pkg_resources.resource_filename('autosolvate', 
    #           os.path.join('data/', solvPrefix, solvent_pdb))
    #   subprocess.call(['cp',solvent_pdb_origin,solvent_pdb])

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
    
    def torchaniFF(self):
        input_pdb = self.solvPrefix + "_solvated.packmol.pdb"
        f_input_pdb = open(input_pdb,'r').readlines()
        line = "%-6s"%"CRYST1" + "%9.3f"%self.cube_size + "%9.3f"%self.cube_size + "%9.3f"%self.cube_size + "%7.2f"%90.0 + "%7.2f"%90.0 + "%7.2f"%90.0 + '%11s'%'P 1'
        print(line)
        f_pdb_pbc = open(self.solvPrefix + "_solvated.packmol.pbc.pdb",'w')
        f_pdb_pbc.write(line+'\n')
        for line in f_input_pdb:
            if 'ATOM' == line.split()[0]:
                f_pdb_pbc.write(line)
        f_pdb_pbc.write('END')
        f_pdb_pbc.close()
        input_pdb_pbc = self.solvPrefix + "_solvated.packmol.pbc.pdb"
       
        atoms = read(input_pdb_pbc)
        print(atoms)
        calculator = torchani.models.ANI1ccx().ase()
        atoms.set_calculator(calculator)
        print("Begin minimizing...")
        opt = BFGS(atoms,trajectory='mmmin.traj')
        opt.run(fmax=0.001)
        print()
        trajmin = Trajectory('mmmin.traj')
        atoms = trajmin[-1]
        atoms.set_calculator(calculator)
        
        ##### mm heat with constant volumn ####
        temp_init = self.tempi
        temp_end = self.temp0
        temp = temp_init
        temp_all = [temp_init]
        
        while temp < temp_end :
         #   print(temp,'K')
            temp = temp + 30
            temp_all.append(temp)
        
        temp_all.append(temp_end)
        Nstep= int(self.nstlim_heating/len(temp_all))
        print(temp_all)
        for temperature in temp_all:
            print('nvt_heating_'+str(temperature)+'.traj')
            dyn = Langevin(atoms, timestep=self.dt * units.fs, temperature_K = temperature , trajectory='nvt_heating_'+str(temperature)+'.traj', logfile='nvt_heating_'+str(temperature)+'.log', friction=0.002)
            dyn.run(Nstep)
            traj = Trajectory('nvt_heating_'+str(temperature)+'.traj')
               # CP = Checkpoint(db='checkpoints_'+temperature+'K.db')
            atoms = traj[-1]
            atoms.set_calculator(calculator)
            traj.close()

                

        

        
        



     #   dyn = NPTBerendsen(atoms, timestep=0.1 * units.fs, temperature_K=300,
      #  taut=10 * units.fs, pressure_au=1.01325 * units.bar,
      #  taup=10 * units.fs, compressibility=4.57e-5 / units.bar,trajectory='npt.traj', logfile='npt.log')
      #  dyn.run(20)
        

        




        


A = MLFF('/home/sunnyxun/projects/MDFF/TorchaniFF/try.xyz')
#A.packSLUSLV()
A.torchaniFF()



