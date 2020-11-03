import pybel
import openbabel as ob
import subprocess

class solventBoxBuilder():
    def __init__(self, xyzfile, solvent='water', slu_netcharge=0, slu_spinmult=1):
        self.xyz = xyzfile
        self.solute = pybel.readfile('xyz', xyzfile).__next__()
        self.slu_netcharge = slu_netcharge
        self.slu_spinmult = slu_spinmult
        # currently hard coded. Can be changed later to be determined automatically based on the density of the solute
        self.solvent = solvent
        self.waterbox_size = 8.0
        # following are for custom organic solvent
        self.cube_size = 26 # in angstrom
        self.slu_pos = self.cube_size/2.0
        self.slv_count = 210
        self.pbcbox_size = self.cube_size + 1

    def getSolutePDB(self):
        print("Converting xyz to pdb")
#        cmd = "babel " + self.xyz + " solute.pdb"
#        subprocess.call(cmd, shell=True)
        obConversion = ob.OBConversion()
        obConversion.SetInAndOutFormats("xyz", "pdb")
        obmol = self.solute.OBMol
        obConversion.WriteFile(obmol,'solute.xyz.pdb')
        # Change the residue name from the default UNL to SLU
        pdb1 = open('solute.xyz.pdb').readlines()
        pdb2 = open('solute.xyz.pdb','w')
        for line in pdb1:
            newline = line.replace('UNL','SLU')
            pdb2.write(newline)
        pdb2.close()
            
    def getFrcmod(self):
        print("Generate frcmod file for the solute")
        print("First generate the gaussian input file")
        cmd1 ="$AMBERHOME/bin/antechamber -i solute.xyz.pdb -fi pdb -o gcrt.com -fo gcrt -gv 1 -ge solute.gesp  -s 2 -nc "+str(self.slu_netcharge) + " -m " + str(self.slu_spinmult)
        print(cmd1)
        subprocess.call(cmd1, shell=True)
        print("Then run Gaussian...")
        cmd2="export GAUSS_SCRDIR=/scratch/$USER/$SLURM_JOBID ;"
        cmd2 = cmd2 + "g16 < gcrt.com > gcrt.out"
        subprocess.call(cmd2, shell=True)
        print("Gaussian ESP calculation done")
        print("Then write out mol2")
        cmd3="$AMBERHOME/bin/antechamber -i solute.gesp -fi gesp -o solute.mol2 -fo mol2 -c resp -eq 2 -rn SLU"
        subprocess.call(cmd3, shell=True)
        print("Finally generate frcmod with parmchk2")
        cmd4 ="$AMBERHOME/bin/parmchk2 -i solute.mol2 -f mol2 -o solute.frcmod"
        subprocess.call(cmd4, shell=True)

    def getHeadTail(self):
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
        print("Now create the solute library file")
        self.writeTleapcmd1()
        cmd ="tleap -s -f leap.cmd > leap_savelib.log"
        subprocess.call(cmd, shell=True)

    def writeTleapcmd_add_water(self):
        print("Now add water  box to the solute")
        self.getHeadTail()
        f = open("leap_add_solventbox.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n")
        f.write("loadamberparams solute.frcmod\n")
        f.write("mol=loadmol2 solute.mol2\n")
        f.write("check mol\n")
        if self.slu_netcharge != 0:
            if self.slu_netcharge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            f.write("addIons2 mol " + ion + " 0\n")
            f.write("check mol\n")
        f.write("solvatebox mol TIP3PBOX " + str(self.waterbox_size) + " iso 0.8  #Solvate the complex with a cubic water box\n")
        f.write("savepdb mol solute_waterbox.pdb\n")
        f.write("saveamberparm mol solute_waterbox.prmtop solute_waterbox.inpcrd\n")
        f.write("quit\n")
        f.close()


    def processPackmolPDB(self):
        data=open('ch3cn_solvated.packmol.pdb').readlines()
        output=open('ch3cn_solvated.processed.pdb','w')
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
        print("Now use packmol to pack the solute and solvent into a box")
        if self.solvent != 'acetonitrile':
            print("The type of solvent is not implemented yet")
            return
        else:
            packmol_inp = open('packmol.inp','w')
            packmol_inp.write("# All atoms from diferent molecules will be at least 2.0 Angstroms apart\n")
            packmol_inp.write("tolerance 2.0\n")
            packmol_inp.write("\n")
            packmol_inp.write("filetype pdb\n")
            packmol_inp.write("\n")
            packmol_inp.write("output ch3cn_solvated.packmol.pdb\n")
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
            packmol_inp.write("structure ch3cn.pdb\n")
            packmol_inp.write("  number " + str(self.slv_count) + " \n")
            packmol_inp.write("  inside cube 0. 0. 0. " + str(self.cube_size) + " \n")
            packmol_inp.write("  resnumbers 2 \n")
            packmol_inp.write("end structure\n")
            packmol_inp.close()

            cmd ="packmol < packmol.inp > packmol.log"
            subprocess.call(cmd, shell=True)
            self.processPackmolPDB()

    def writeTleapcmd_ch3cn_solvated(self):
        f = open("leap_packmol_solvated.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
        f.write("loadamberparams ch3cn.frcmod\n")
        f.write("loadAmberPrep ch3cn.prep\n")
        f.write("loadamberparams solute.frcmod\n")
        f.write("loadoff solute.lib\n")
        f.write("\n")
        f.write("SYS = loadpdb ch3cn_solvated.processed.pdb\n")
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
        f.write("saveamberparm SYS ch3cn_solvated.prmtop ch3cn_solvated.inpcrd #Save AMBER topology and coordinate files\n")
        f.write("savepdb SYS ch3cn_solvated.pdb\n")
        f.write("quit\n")
        f.close()

    def createAmberParm(self):
        print("Generate Amber parameters for the solvated system")
        if self.solvent == 'water':
            self.writeTleapcmd_add_water()
            cmd ="tleap -s -f leap_add_solventbox.cmd > leap_add_solventbox.log"
            subprocess.call(cmd, shell=True)
        elif self.solvent == 'acetonitrile':
            self.packSLUSLV()
            self.writeTleapcmd_ch3cn_solvated()
            cmd ="tleap -s -f leap_packmol_solvated.cmd > leap_packmol_solvated.log"
            subprocess.call(cmd, shell=True)

if __name__ == '__main__':
    builder = solventBoxBuilder('coumarin7oh_min.xyz', 'acetonitrile', slu_netcharge=1, slu_spinmult=2)
    builder.getSolutePDB()
    builder.getFrcmod()
    builder.createLib()
    builder.createAmberParm()
