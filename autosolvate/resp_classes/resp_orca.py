import autosolvate.resp_classes.gen_esp as gen_esp
from autosolvate.resp_classes.resp_abstract import RespABC
import os
import subprocess
#import os, re, subprocess, shutil, glob, sys
def get_crds_from_xyz(inp):
    crds = []
    with open(inp,'r') as f:
        data = f.readlines()
        for line in data[2:]:
            if len(line) > 3:
                line = line.strip()
                atom = line.split()[0].capitalize()
                x = line.split()[1]
                y = line.split()[2]
                z = line.split()[3]
                crds.append([atom,x,y,z])
    return crds


class RespORCA(RespABC):
    def __init__(self, **kwargs):
        print("*"*40)
        print("Run Orca to generate RESP charge and mol2 file".center(40," ") )
        print("*"*40)
        self.srun_use = kwargs["srun_use"] if "srun_use" in kwargs else False
        # TODO: read in srun_use option from kwargs
        RespABC.__init__(self, **kwargs)
        
        if self.qm_exe == None:
           print("WARNING: Orca executable name is not specified for RESP charge fitting!")
           print("WARNING: Using orca by default. If failed later,\
                  please rerun with the option -e specified!")
           self.qm_exe = 'orca'
        if self.qm_dir == None:
            print("WARNING: orca executable directory is not specified for RESP charge fitting!")
            print("WARNING: Setting to default path: /opt/orca/5.0.2/")
            print("WARNING: If failed later, please rerun with the option -d specified!")
            self.qm_dir = '/opt/orca/5.0.2/'
        
        self.orcapath = os.path.join(self.qm_dir, self.qm_exe)

    def writeOrcaInput(self):
        """
        Set up Orca calculation to compute electrostatic potential.
     
        optional arguments:
          calculation - one of 'optimize' (hours), or 'singlepoint' (minutes) (defaults to 'singlepoint')
     
        """
        print("-"*40)
        print("Preparing Orca input files...".center(40," "))
        print("-"*40)

        charge = self.molecule.GetTotalCharge()
        multiplicity = self.molecule.GetTotalSpinMultiplicity()

        orca_inp = self.molname + "_orca.in"
        orca_out = self.molname + '_orca.out'
        xyzfile = self.xyzfile
        orca_crd = get_crds_from_xyz(xyzfile)
        
        with open(orca_inp,'w') as ofile:
            ofile.write('! b3lyp 6-31g* keepdens\n\n')
            ofile.write("%pal\nnprocs " + str(self.nprocs) + "\nend\n\n")
            ofile.write('* xyz ' + str(charge) + ' ' + str(multiplicity) + '\n')
            for line in orca_crd:
                atomname = line[0]
                x  = line[1]
                y = line[2]
                z = line[3]
                ofile.write(atomname + '  ' + x + '  ' + y  + '  ' + z + '\n')
            ofile.write('*')
            
        cmd1 = self.orcapath + ' ' + orca_inp + ' > ' + orca_out
        
        ## check convergence of orcaout ##
        convergence = False
        if os.path.exists(orca_out):
            with open(orca_out) as f:
                data = f.read()
                if 'ORCA TERMINATED NORMALLY' in data:
                    print('Orca calculation is finished for esp...')
                    convergence = True
            
        if convergence == False:
            if self.srun_use:
                cmd1='srun -n ' + self.nprocs + ' ' +cmd1
            print(cmd1)
            subprocess.call(cmd1, shell=True)
    
    def generateESP(self):
        vpotpath = self.orcapath + '_vpot'
        orca_inp = self.molname + "_orca.in"
        orca_out = self.molname + '_orca.out'
        gbw = self.molname + '_orca.gbw'
        density = self.molname + '_orca.scfp'
        espf =  self.molname + '_orca.esp'
        out =  self.molname + '_orca.espout'
        crds = get_crds_from_xyz(self.xyzfile)
        elements = []
        newcrds = []
        for line in crds:
            newcrds.append([float(line[1]),float(line[2]),float(line[3])])
            elements.append(line[0].capitalize())
        
        gen_esp.gen_grids(crds=newcrds,elements=elements,orcapath=vpotpath,gbw=gbw,denisty=density,out=out)
        gen_esp.convertoesp(espin=out,espout=espf,crds=newcrds)
        
       # os.system("resp -O -i resp1.in -o resp1.out -p resp1.pch -t resp1.chg \
       #        -e %s -s resp1_calc.esp" %espf)
      #  os.system("resp -O -i resp2.in -o resp2.out -p resp2.pch -q resp1.chg \
       #       -t resp2.chg -e %s -s resp2_calc.esp" %espf)
        
    
    
    def respFit(self):
        Method1 = True
        print("-"*40)
        print("Start RESP fitting with ESP generaged by orca".center(40," ") )
        print("-"*40)
        
        ac_filename = self.molname + ".ac"
        print("Generating Antechamber file: " + ac_filename)
        
        cmd = "antechamber -i {pdbfile} -fi pdb -o {ac} -fo ac -c dc -nc {charge} -pl 30".format(pdbfile=self.pdbfile, ac=ac_filename, charge=self.charge)
        print(cmd)
        with open('antechamber_acfile.log', 'w') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                           
        resp1_filename = self.molname + '.resp1'
        print("Using respgen to generate stage 1 RESP fitting input file: " + resp1_filename)
        
        cmd = "respgen -i {ac} -o {resp1} -f resp1 -l 10".format(ac=ac_filename, resp1=resp1_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)

        # Stage 2
        resp2_filename = self.molname + '.resp2'
        print("Using respgen to generate stage 2 RESP fitting input file: " + resp2_filename)
        
        cmd = "respgen -i {ac} -o {resp2} -f resp2 -l 10".format(ac=ac_filename, resp2=resp2_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)
        
        esp_filename =  self.molname + '_orca.esp'
        
        # Perform RESP fit.
        #Stage 1
        qout1_filename = self.molname + '.qout_stage1'
        respout1_filename = self.molname + '.respout1'
        print("Using resp to perform stage 1 RESP fitting to generate: " + qout1_filename)
        
        cmd = "resp -O -i {resp1} -o {respout1} -e {esp} -t {qout1}".format(
            resp1=resp1_filename, respout1=respout1_filename,
            esp=esp_filename, qout1=qout1_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)

        #Stage 2
        qout2_filename = self.molname + '.qout_stage2'
        respout2_filename = self.molname + '.respout2'
        print("Using resp to perform stage 2 RESP fitting to generate: " + qout2_filename)

        cmd = "resp -O -i {resp2} -o {respout2} -e {esp} -q {qout1} -t {qout2}".format(
            resp2=resp2_filename, respout2=respout2_filename, esp=esp_filename,
            qout1=qout1_filename, qout2=qout2_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)

        # Write molecule with updated charges to mol2 file.
        mol2_filename = self.molname + ".mol2"
        print("Writing out the mol2 file with resp charge: " + mol2_filename)
        cmd = "antechamber -i {ac} -fi ac -o {mol2} -fo mol2 -c rc -cf {qout2} -pl 30".format(
            ac=ac_filename, mol2=mol2_filename, qout2=qout2_filename)
        print(cmd)
        subprocess.call(cmd,shell=True)
        
 
         #   with open('resp1.in','w') as fresp1:
          #      fresp1.write()
            
    def run(self):
      #  if not os.path.isdir(self.resp_scr_dir):
       #     print("Creating the scratch folder for RESP fitting: ", self.resp_scr_dir)
       #     os.mkdir(self.resp_scr_dir)

      #  print("Copying over the pdb file", self.pdbfile, " to ", self.resp_scr_dir)
      #  shutil.copy(os.path.join(self.rundir,self.pdbfile), self.resp_scr_dir)

      #  print("Navigating to the scratch folder for RESP fitting: ", self.resp_scr_dir)
      #  os.chdir(self.resp_scr_dir)
        
        self.writeOrcaInput()
        self.generateESP()
        self.respFit()

     #   print("Navigating back to the folder for AutoSolvate run: ", self.rundir)
     #   os.chdir(self.rundir)
        
     #   mol2_filename = self.molname + ".mol2"
    #    print("Copying the generated mol2 file ", mol2_filename, " to rundir", self.rundir)
     #   shutil.copy(os.path.join(self.resp_scr_dir, mol2_filename), self.rundir)

        
        
        
        
        
        
        
        
        

        
        
        