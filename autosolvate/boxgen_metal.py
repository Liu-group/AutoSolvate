#import autosolvate
import getopt, sys, os
import subprocess
from glob import glob
import pkg_resources

Number_A = 6.02214076 * 10**23

amber_solv_dict = {'water':     [' ','TIP3PBOX '],
                   'methanol':  ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
                   'chloroform':['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
                   'nma':       ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}

custom_solv_dict = {'acetonitrile':'ch3cn'}
custom_solv_residue_name = {'acetonitrile':'C3N'}
elements_dic = {'H':1,'HE':2,'LI':3,'BE': 4,'B':5,'C':6, 'N':7, 'O':8, 'F':9,'NE':10,
                 'NA': 10,'MG':12,'AL':13,'SI':14,'P':15,'S':16, 'CL':17, 'AR':18,'K':19,'CA':20}

weight_dic ={'water': 18.01528 ,'acetonitrile': 41.05, 'methanol':32.04, 'chloroform':119.38, 'nma':107.156 } ##g/mol

density_dic = {'water' : 997,'acetonitrile': 786, 'methanol':792, 'chloroform': 1490,'nma':990} ##kg/m3

closeness_dic = {'water': 0.5, 'acetonitrile':1.80, 'methanol':0.6, 'nma': 0.58 ,'chloroform':0.58}

def calculateSolvnumber(solvPrefix,volume):  ###Volume should be m3
    denisty = density_dic[solvPrefix]
    weight = weight_dic[solvPrefix]
    mass = volume * denisty  ## kg
    Mol = mass*1000/weight
    NumberOfSolv = Mol * Number_A
    return int(NumberOfSolv) 

class solventBoxBuilderMetal(object):
    def __init__(self, pdb_prefix, totalcharge, 
                 solvent , solvent_frcmod, solvent_off,
                 slv_count, 
                 cube_size , closeness,  outputFile, amberhome):
        self.pdb_prefix = pdb_prefix
        self.solvent = solvent
        self.totalcharge = totalcharge
        self.solvent_frcmod=solvent_frcmod
        self.solvent_off = solvent_off
        self.cube_size = cube_size
        self.slu_pos = self.cube_size/2.0
        self.pbcbox_size = self.cube_size+2
        self.amberhome = amberhome
        self.closeness = closeness
        

        if closeness in [ 'automated','default','Default']:
            if self.solvent in ['acetonitrile','ch3cn','']:
                self.closeness=1.88
            elif self.solvent=='water':
                self.closeness=0.50
            elif self.solvent=='methanol':
                self.closeness=0.60
            elif self.solvent=='nma':
                self.closeness=0.58
            elif self.solvent=='chloroform':
                self.closeness=0.58
            else:
                print("Warning: solvent not supported for automated closeness determination")
        else:
            self.closeness = 0.5
        self.waterbox_size = 8.0
        self.cube_size = cube_size
        self.volumne = (self.cube_size * 10**(-10))**3
        if self.solvent not in ['acetonitrile','water','methanol','nma','chloroform']:
            print("Warning: solvent not supported for automated counting determination, using the default number")
            slv_count = 210*8
            self.slv_count = slv_count
        else:
            self.slv_count = calculateSolvnumber(self.solvent, self.volumne)
        
        self.outputFile=outputFile
        if self.outputFile == "":
            self.outputFile = self.pdb_prefix + "_solvated"
        
        if self.totalcharge.upper() in ['DEFAULT']:
            if self.pdb_prefix + '.info' in glob('*.info'):
                with open(self.pdb_prefix + '.info', 'r') as f:
                    data = f.readlines()

                begin, end = 0, len(data) 
                for sort, line in enumerate(data):
                    if '@charges' in line:
                        begin = sort + 1
                    elif '@connection_infos' in line:
                        end = sort
                newtotalcharge = 0  

                for line in data[begin:end]:
                    if 'LG' not in line:
                        parts = line.split()
                        if len(parts) == 2:
                            metalcharge = int(parts[1])
                            newtotalcharge += metalcharge
                    else:
                        parts = line.split()
                        if len(parts) == 2:
                            ligandchg = int(parts[1])
                            newtotalcharge += ligandchg
                self.totalcharge = newtotalcharge
            else:
                print('Error: can not find',self.pdb_prefix + '.info', 'please assign the charge manually') 
                sys.exit()
        else:
            try:
                self.totalcharge = int(totalcharge)
               # print('The totalcharge is assigned as',self.totalcharge)
            except ValueError:
                print('Error: Invalid input for the number of totalcharge')
                sys.exit()         

    def get_ligand_name(self):
        lignames = []
        pdbfile = self.pdb_prefix + '_mcpbpy.pdb'
        if pdbfile not in glob("*.pdb"):
            print('Error:',pdbfile,'can not be found, please check the previous steps!')
            sys.exit()
        with open (pdbfile,'r') as f:
            for line in f:
                if len(line) > 20:
                    if line.split()[0] in ['ATOM',"HETATM"]:
                        resname = line[17:20].strip()
                        if resname not in lignames:
                            lignames.append(resname)
        return lignames
    
    def get_frcmond_name(self):
        frcmondnames = []
        for i in glob('LG*.frcmod'):
            frcmondnames.append(i)
        return frcmondnames
    
    def writeTleapcmd_add_solvent(self,lignames,frcmondnames):
        if self.pdb_prefix + '_tleap.in' in glob('*'):
            tleapin = open(self.pdb_prefix + '_tleap.in', 'r')
          #  print(self.pdb_prefix + '_tleap.in')
            data_tleapin = tleapin.readlines()
            tleapin.close()
            ofile = open("leap_add_solventbox.cmd","w")
         #   print('leap_add_solventbox.cmd')
            ofile.write("source leaprc.protein.ff14SB\n")
            ofile.write("source leaprc.gaff\n")
            ofile.write("source leaprc.water.tip3p\n")
            ofile.write(str(amber_solv_dict[str(self.solvent)][0]))
            ofile.write('addAtomTypes {\n')
            AtomTypes = []
            bondmol = []
            for line in data_tleapin:
                if len(line) > 1:
                    if '{' in line:
                        if '}' in line:
                            if line not in AtomTypes:
                                AtomTypes.append(line) 
                if "bond mol" in line:
                    if line not in bondmol:
                        bondmol.append(line)
            
            for i in AtomTypes:
                ofile.write(i)
            
            ofile.write('}\n')
            
            for lg in lignames:
                ofile.write(lg + ' = loadmol2 ' + lg + '.mol2\n')
            
            for frcmod in frcmondnames:
                ofile.write('loadamberparams ' + frcmod + '\n')
            ofile.write('loadamberparams frcmod.ionslm_126_opc\n')
            ofile.write('loadamberparams ' + self.pdb_prefix + '_mcpbpy.frcmod\n')
            ofile.write('mol = loadpdb ' + self.pdb_prefix + '_mcpbpy.pdb\n')

            for i in bondmol:
             #   print('bond mol',i )
                ofile.write(i)
            
            ofile.write("solvatebox mol " + (str(amber_solv_dict[str(self.solvent)][1])) + str(self.slu_pos) + " iso "+str(self.closeness)+"  #Solvate the complex with a cubic solvent box\n")
            if self.totalcharge != 0:
                if self.totalcharge > 0:
                    ion = 'Cl-'
                else:
                    ion = 'Na+'
                
                ionnum = str(abs(self.totalcharge))
                ofile.write("addIons2 mol " + ion + ' ' + ionnum + '\n')
            
            ofile.write("check mol\n")
            ofile.write("savepdb mol " + str(self.outputFile) + ".pdb\n")
            ofile.write("saveamberparm mol " + str(self.outputFile) + ".prmtop " + str(self.outputFile) + ".inpcrd\n")
            ofile.write("quit\n")
            ofile.close()         
        
        else:
            print('Error: can not find',self.pdb_prefix + '_tleap.in','Please check the files in previous steps')
            sys.exit()
    
    def processPackmolPDB(self,output_pdb):
        r"""
        Convert file for custom solvents like CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        data=open(output_pdb,'r').readlines()
        output=open(self.pdb_prefix + '_processed.pdb','w')
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
        solutePDB = self.pdb_prefix + '_mcpbpy.pdb'
        print("Now use packmol to pack the solute and solvent into a box")
        if self.solvent not in  custom_solv_dict.keys():
            print("The type of solvent is not implemented yet")
            sys.exit()
        else:
            solvPrefix = custom_solv_dict[self.solvent]
            solvent_pdb = solvPrefix+'.pdb'
            solvent_pdb_origin = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data/', solvPrefix, solvent_pdb))
            subprocess.call(['cp',solvent_pdb_origin,solvent_pdb])
      #      init_file_path = autosolvate.__file__
      #      autosolvate_dir = os.path.dirname
      #      solvent_pdb_origin = os.path.join(autosolvate_dir, 'data', solvPrefix, solvent_pdb)
      #      print(solvent_pdb_origin)
     #       subprocess.call(['cp',solvent_pdb_origin,solvent_pdb])
            output_pdb = self.pdb_prefix + "_packmol.pdb"
            packmol_inp = open('packmol.inp','w')
            packmol_inp.write("# All atoms from diferent molecules will be at least %s Angstroms apart\n" % self.closeness)
            packmol_inp.write("tolerance %s\n" % self.closeness)
            packmol_inp.write("\n")
            packmol_inp.write("filetype pdb\n")
            packmol_inp.write("\n")
            packmol_inp.write("output " + output_pdb + "\n")
            packmol_inp.write("\n")
            packmol_inp.write("# add the solute\n")
            packmol_inp.write("structure " + solutePDB +  "\n")
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
            print(output_pdb)
            self.processPackmolPDB(output_pdb)

    def writeTleapcmd_add_solvent_custom(self,lignames,frcmondnames):
        if self.pdb_prefix + '_tleap.in' in glob('*'):
            tleapin = open(self.pdb_prefix + '_tleap.in', 'r')
            data_tleapin = tleapin.readlines()
            tleapin.close()

            inputpdb = self.pdb_prefix + '_processed.pdb'
            solvPrefix = custom_solv_dict[self.solvent]
            solvent_pdb = solvPrefix+'.pdb'

          #  autosolvate_dir = os.path.dirname(autosolvate.__file__)
         #   solvent_frcmod = solvPrefix + '.frcmod'
          #  solvent_frcmod_path = os.path.join(autosolvate_dir, 'data', solvPrefix, solvent_frcmod)

          #  solvent_prep = solvPrefix + '.prep'
          #  solvent_prep_path = os.path.join(autosolvate_dir, 'data', solvPrefix, solvent_prep)
       
            solvent_frcmod = solvPrefix+'.frcmod'
            solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data',solvPrefix,solvent_frcmod))

            solvent_prep = solvPrefix+'.prep'
            solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data',solvPrefix,solvent_prep))
         
            f = open("leap_add_solventbox.cmd","w")
            f.write("source leaprc.protein.ff14SB\n")
            f.write("source leaprc.gaff\n")
            f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
            f.write('addAtomTypes {\n')
            AtomTypes = []
            bondmol = []
            for line in data_tleapin:
                if len(line) > 1:
                    if '{' in line:
                        if '}' in line:
                            if line not in AtomTypes:
                                AtomTypes.append(line) 
                if "bond mol" in line:
                    if line not in bondmol:
                        bondmol.append(line)
            
            for i in AtomTypes:
                f.write(i)
            
            f.write('}\n')
            
            for lg in lignames:
               f.write(lg + ' = loadmol2 ' + lg + '.mol2\n')

            
            for frcmod in frcmondnames:
                f.write('loadamberparams ' + frcmod + '\n')
            f.write('loadamberparams frcmod.ionslm_126_opc\n')
            f.write("loadamberparams " + solvent_frcmod_path + "\n")
            f.write("loadAmberPrep " + solvent_prep_path + "\n")
            f.write("loadamberparams "+ self.pdb_prefix + "_mcpbpy.frcmod"+"\n")
            f.write("mol = loadpdb " + inputpdb +" \n")
            for i in bondmol:
                f.write(i)
            f.write("\n")
            
            f.write("check mol\n")
            f.write("\n")
            if self.totalcharge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            if self.totalcharge != 0:
                if self.totalcharge > 0:
                    ion = 'Cl-'
                else:
                    ion = 'Na+'
                f.write("addIons2 mol " + ion + ' ' + str(abs(self.totalcharge)) +" \n")
                f.write("check mol\n")
            f.write("# set the dimension of the periodic box\n")
            f.write("set mol box {" + str(self.pbcbox_size) +", " + str(self.pbcbox_size) + ", " + str(self.pbcbox_size) + "}\n")
            f.write("\n")
            f.write("saveamberparm mol "+str(self.outputFile)+".prmtop "+str(self.outputFile)+".inpcrd #Save AMBER topology and coordinate files\n")
            f.write("savepdb mol "+str(self.outputFile)+".pdb\n")
            f.write("quit\n")
            f.close()
        else:
            print('Error: can not find',self.pdb_prefix + '_tleap.in','Please check the files in previous steps')
            sys.exit()

    def writeTleapcmd_add_solvent_by_user(self,lignames,frcmondnames):
        if self.pdb_prefix + '_tleap.in' in glob('*'):
            tleapin = open(self.pdb_prefix + '_tleap.in', 'r')
            data_tleapin = tleapin.readlines()
            tleapin.close()
            ofile = open("leap_add_solventbox.cmd","w")
            ofile.write("source leaprc.protein.ff14SB\n")
            ofile.write("source leaprc.gaff\n")
            ofile.write("source leaprc.water.tip3p\n")
            ofile.write('addAtomTypes {\n')
            AtomTypes = []
            bondmol = []
            for line in data_tleapin:
             #   print(line)
                if len(line) > 1:
                    if '{' in line:
                        if '}' in line:
                            if line not in AtomTypes:
                                AtomTypes.append(line) 
                if "bond mol" in line:
                    if line not in bondmol:
                        bondmol.append(line)
            
            for i in AtomTypes:
                ofile.write(i)
            
            ofile.write('}\n')

            
            for lg in lignames:
                ofile.write(lg + ' = loadmol2 ' + lg + '.mol2\n')
            
            for frcmod in frcmondnames:
                ofile.write('loadamberparams ' + frcmod + '\n')
        
                print('loadamberparams ' + frcmod + '\n')
            ofile.write('loadamberparams frcmod.ionslm_126_opc\n')
            ofile.write('loadamberparams ' + self.pdb_prefix + '_mcpbpy.frcmod\n')
            ofile.write('mol = loadpdb ' + self.pdb_prefix + '_mcpbpy.pdb\n')

            for i in bondmol:
                ofile.write(i)

            print("Now add custom pre-equlibrated solvent box to the solute")
            ofile.write("loadoff " + self.solvent_off + "\n")
            ofile.write("loadamberparams " + self.solvent_frcmod + "\n")
            ofile.write("loadamberparams "+ self.pdb_prefix + "_mcpbpy.frcmod"+"\n")
            
            ofile.write("solvatebox mol " + self.solvent + " " 
                    + str(self.slu_pos) + " iso 0.8  #Solvate the complex with a cubic solvent box\n") 
            # Notice that we want to add the ion after solvation because we don't want the counter ion to be too close to solute
            if self.totalcharge != 0:
                if self.totalcharge > 0:
                    ion = 'Cl-'
                else:
                    ion = 'Na+'
                ofile.write("addIons2 mol " + ion + " 0\n")
                ofile.write("check mol\n")
            ofile.write("check mol\n")
            ofile.write("savepdb mol " + str(self.outputFile) + ".pdb\n")
            ofile.write("saveamberparm mol " + str(self.outputFile) + ".prmtop " + str(self.outputFile) + ".inpcrd\n")
            ofile.write("quit\n")
            ofile.close()
    
    def tleap(self):
        cmd = self.amberhome + 'tleap -f ' + 'leap_add_solventbox.cmd > tleap.log'
        subprocess.call(cmd,shell=True)
        
    def build(self):
        frcmondnames = self.get_frcmond_name()
        lgnames = self.get_ligand_name()
        if self.solvent in custom_solv_dict:
            self.packSLUSLV()
            self.writeTleapcmd_add_solvent_custom(lignames=lgnames,frcmondnames=frcmondnames)
        elif len(self.solvent_frcmod) > 0 and len(self.solvent_off) > 0:
            self.writeTleapcmd_add_solvent_by_user(lignames=lgnames,frcmondnames=frcmondnames)
        if self.solvent in amber_solv_dict.keys():
            self.writeTleapcmd_add_solvent(lignames=lgnames,frcmondnames=frcmondnames)
        self.tleap()

def startboxgen(argumentList):
    options = "hm:s:o:c:b:a:t:l:p:"
    long_options = ["help", "pdb_prefix", "solvent", "output", "totalcharge", 
                "cubesize", "amberhome","closeness","solventoff","solventfrcmod"]
    arguments, values = getopt.getopt(argumentList, options, long_options)

    pdb_prefix= '' 
    totalcharge='Default' 
    solvent = "water"
    solvent_frcmod = ""
    solvent_off = ""
    slv_count = 210*8
    cube_size = 54
    closeness = "automated"
    outputFile = ""
    amberhome = '$AMBERHOME/bin/'
    for currentArgument, currentValue in arguments:
        if  currentArgument in ("-h", "--help"):
            print('Usage: autosolvate_metal boxgen [OPTIONS]')
            print('  -m, --pdb_prefix           prefix of pdb file name not include _mcpb')
            print('  -s, --solvent              name of solvent')
            print('  -o, --output               prefix of the output file names')
            print('  -c, --charge               formal charge of solute')
            print('  -b, --cubesize             size of solvent cube in angstroms')
            print('  -a, --amberhome            path to the AMBER molecular dynamics package root directory')
            print('  -t, --closeness            Solute-solvent closeness setting')
            print('  -l, --solventoff           path to the custom solvent .off library file')
            print('  -p, --solventfrcmod        path to the custom solvent .frcmod file')
            print('  -h, --help                 short usage description')
            exit()
        elif currentArgument in ('-m','-pdb_prefix'):
            pdb_prefix = str(currentValue)
        elif currentArgument in ("-s", "--solvent"):
            solvent=str(currentValue)
        elif currentArgument in ("-o", "--output"):
            print ("Output:", currentValue)
            outputFile=str(currentValue)
        elif currentArgument in ("-c", "--charge"):
            print ("Charge:", currentValue)
            totalcharge=str(currentValue)
        elif currentArgument in ("-b", "--cubesize"):
            print ("Cubesize:", currentValue)
            cube_size=float(currentValue)
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
    builder = solventBoxBuilderMetal(pdb_prefix=pdb_prefix,solvent=solvent,
                                        outputFile=outputFile,totalcharge=totalcharge,
                                        cube_size=cube_size,amberhome=amberhome,closeness=closeness,
                                        solvent_frcmod=solvent_frcmod,slv_count=slv_count,solvent_off=solvent_off)
    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startboxgen(argumentList)



     
    

    



        