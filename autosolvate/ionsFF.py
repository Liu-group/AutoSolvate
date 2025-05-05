from pymsmt.mol.element import get_ionljparadict
from pymsmt.mol.element import METAL_PDB
import sys
import subprocess
import os
import re
import pkg_resources

Number_A = 6.02214076 * 10**23

OrganicMass = { 'H' : 1.008,
                'C' : 12.01,
                'N' : 14.01,
                'O' : 16.00,
                'S' : 32.06,
                'P' : 30.97,
              }
organic_LJ  = {'S':{'sigma':2.100, 'epsilon':0.2500},'P':{'sigma':2.100, 'epsilon':0.2000},'N':{'sigma':1.8240, 'epsilon':0.1700}} #### from GAFF force field

amber_solv_dict = {'water': [' ','TIP3PBOX '],
                   'methanol': ['loadOff solvents.lib\n loadamberparams frcmod.meoh\n', 'MEOHBOX '],
                   'chloroform': ['loadOff solvents.lib\n loadamberparams frcmod.chcl3\n', 'CHCL3BOX '],
                   'nma': ['loadOff solvents.lib\n loadamberparams frcmod.nma\n', 'NMABOX ']}

custom_solv_dict = {'acetonitrile':'ch3cn'}
custom_solv_residue_name = {'acetonitrile':'C3N'}

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


def write_filtered_mol2(input_file, output_file, custom_line):
    with open(input_file, "r") as f:
        lines = f.readlines()
    
    skip_next = False  # 标记是否跳过下一行
    filtered_lines = []  # 存储要写入的行

    for line in lines:
        if skip_next:
            # 在跳过原本的下一行后，写入自定义的新行
            filtered_lines.append(custom_line + "\n")
            skip_next = False  # 恢复标记
            continue  # 跳过被替换的行

        if line.strip() == "@<TRIPOS>ATOM":
            skip_next = True  # 发现 @<TRIPOS>ATOM，下一行要被替换

        filtered_lines.append(line)  # 其他行正常存入

    # 写入新文件
    with open(output_file, "w") as f:
        f.writelines(filtered_lines)


class file_prep_for_ion():
    def __init__(self, xyzfile, charge,solvent,cubic_size,closeness,solvent_off="", solvent_frcmod=""):
        self.closeness = closeness
        self.cube_size = cubic_size
        self.xyzfile = xyzfile
        self.charge = charge
        self.solvent = solvent
        self.boxsize = cubic_size
        self.slu_pos = str(float(self.cube_size)/2.0)
        self.volumne = (float(self.cube_size) * 10**(-10))**3
        self.solvent_off = solvent_off
        self.solvent_frcmod = solvent_frcmod
        
        if self.solvent not in ['acetonitrile','water','methanol','nma','chloroform']:
            print("Warning: solvent not supported for automated counting determination, using the default number")
            slv_count = 210*8
            self.slv_count = slv_count
        else:
            self.slv_count = calculateSolvnumber(self.solvent, self.volumne)
        
       # print(self.xyzfile)
        
    def get_mass(self):
        with open(self.xyzfile, 'r') as f:
            lines = f.readlines()
            num_atoms = int(lines[0].strip())
            print('Number of atoms in the system:', num_atoms)

            if num_atoms != 1:
                raise ValueError("This function only supports single-atom systems.")

            print('This is a single atom system')
            self.element_symbol = lines[2].split()[0].strip().upper()
            print('Element symbol:', self.element_symbol)

            # --- Step 1: 优先匹配 key[0] ---
            matched_key = next((key for key in METAL_PDB if key[0] == self.element_symbol), None)

            # --- Step 2: 如果 key[0] 匹配不到，再尝试 key[1] ---
            if not matched_key:
                matched_key = next((key for key in METAL_PDB if key[1] == self.element_symbol), None)

            if matched_key:
                print('This is a metal ion!')
                print(f'Matched key: {matched_key} ')

                try:
                    atomnumber = METAL_PDB[matched_key][0]
                    mass = METAL_PDB[matched_key][1]
                    self.atomnumber = atomnumber
                    self.mass = mass
                    self.metal_type = matched_key[1]  # 可选：记录真实金属类型（如 CU, FE2）
                    print(f"Mass = {self.mass}, Atom Number = {self.atomnumber}, Type = {self.metal_type}")
                    return 'metal'
                except Exception as e:
                    raise ValueError(f"Error retrieving data for metal ion {matched_key}: {e}")

            # --- Step 3: 检查是否是有机离子 ---
            print(f"{self.element_symbol} is not a supported ion in MCPB.")
            print("Try to search organic ions...")

            if self.element_symbol in ['H', 'C', 'N', 'O', 'P', 'S']:
                print("This is an organic ion.")
                try:
                    self.mass = OrganicMass[self.element_symbol]
                    return 'organic'
                except KeyError:
                    raise ValueError(f"OrganicMass dictionary does not contain {self.element_symbol}")
            else:
                raise TypeError(
                    f"Error: This is not a metal ion or organic ion supported by the current version of Autosolvate.\n"
                    f"Element symbol: {self.element_symbol}"
            )


        
    def gen_LJ_for_ions(self):
        IonLJParaDict = get_ionljparadict('opc')
        keyname = self.element_symbol.capitalize() + str(int(self.charge))
        
        if keyname in IonLJParaDict.keys():
            print('This is a supported ion in MCPB')
            self.LJ = IonLJParaDict[keyname]
            self.sigma = self.LJ[0]
            self.epsilon = self.LJ[1]
        else:
            if len(keyname) == 2:
                newkeyname = keyname + ' '
                if newkeyname in IonLJParaDict.keys():
                    self.LJ = IonLJParaDict[newkeyname]
                    self.sigma = self.LJ[0]
                    self.epsilon = self.LJ[1]
                else:                
                    print('This is not a supported ion in MCPB')
                    print('try to search organic ions')
                    if self.element_symbol in [ 'S', 'P', 'N']:
                        print('This is an organic ion')
                        self.LJ = organic_LJ[self.element_symbol]
                        self.sigma = self.LJ['sigma']
                        self.epsilon = self.LJ['epsilon']

                    else:
                        raise TypeError('Error: This is not a metal ion or organic ion supported by current version of Autosolvate, please check the input file and charge')
    
    def write_mol2(self):
        pdb = os.path.splitext(self.xyzfile)[0] + '.pdb'
        cmd = 'obabel '+ self.xyzfile + ' -O '+ pdb
        subprocess.call(cmd,shell=True)
        
        with open(pdb,'r') as f:
            lines = f.readlines()
            for line in lines:
                if line.split()[0] in ['ATOM','HETATM']:
                    resname = line[17:20].strip()
                    atomname = line[12:16].strip()
                    self.resname = resname
                    self.atomname = atomname
                    break
        
        mol2 = os.path.splitext(self.xyzfile)[0] + '.mol2'
        cmd = '$AMBERHOME/bin/metalpdb2mol2.py -i ' + pdb + ' -o ' +  mol2 + ' -c ' + str(self.charge)
        subprocess.call(cmd,shell=True)
        self.mol2 = mol2
    
    def write_frcmod(self):
        frcmod = os.path.splitext(self.xyzfile)[0] + '.frcmod'
        with open(frcmod,'w') as f:
            f.write('MASS\n')
            f.write('  '+self.atomname+'  '+str(self.mass)+'\n')
            f.write('NONB\n')
            f.write('  '+self.atomname+'  '+str(self.sigma)+'  '+str(self.epsilon)+'\n')
    
    def process_mol2(self):
        processed_mol2 = os.path.splitext(self.xyzfile)[0] + '_processed.mol2'
        with open(self.mol2,'r') as f:
            lines = f.readlines()
            for count, line in enumerate(lines):
                if line.strip() == '@<TRIPOS>ATOM':
                    atomline = lines[count+1]
                    pattern = self.atomname
                    re_pattern = rf'(\s){pattern}(\s+)'
                    
                    def replace_second_match():
                        matchline = 0  # 计数器
                        def inner_replace(match):
                            nonlocal matchline  
                            matchline += 1
                            if matchline == 2:  # 只替换第二个匹配到的
                                return f"{match.group(1)}X1{match.group(2)}"
                            return match.group(0)  # 保持不变
                        return inner_replace

                    updated_line = re.sub(re_pattern, replace_second_match(), atomline)
        
        write_filtered_mol2(self.mol2, processed_mol2, updated_line)
        self.process_mol2 = processed_mol2
    
    def write_frcmod(self):
        frcmod = os.path.splitext(self.xyzfile)[0] + '.frcmod'
        self.frcmod = frcmod
        with open(frcmod,'w') as f:
            f.write('REMARK GOES HERE, THIS FILE IS GENERATED BY AUTOSOLVATE\n')
            f.write('MASS\n')
            f.write('  X1  '+str(self.mass)+'\n\n')
            f.write('BOND\n\n')
            f.write('ANGLE\n\n')
            f.write('DIHE\n\n')
            f.write('IMPR\n\n')
            f.write('NONB\n')
            f.write('  X1          '+str(self.sigma)+'  '+str(self.epsilon)+'\n')
    
    def writeTleapcmd_add_solvent(self):
        r"""
        Write tleap input file to add solvent

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        ofile = open("leap_add_solventbox.cmd","w")
        ofile.write("source leaprc.protein.ff14SB\n")
        ofile.write("source leaprc.gaff\n")
        ofile.write("source leaprc.water.tip3p\n")
        ofile.write(str(amber_solv_dict[str(self.solvent)][0]))
        symbol = 'X1'
        hybrid = 'sp3'
        output = f'addAtomTypes {{\n    {{ "{symbol}" "{self.atomname}" "{hybrid}" }}\n}}'
        ofile.write(output)
        ofile.write('\n')
        ofile.write(self.resname + '= loadmol2 '+self.process_mol2+'\n') 
        ofile.write('loadamberparams '+self.frcmod+'\n')
        ofile.write('mol = loadmol2 '+os.path.splitext(self.xyzfile)[0] + '_processed.mol2\n')
        ofile.write("solvatebox mol " + (str(amber_solv_dict[str(self.solvent)][1])) + str(self.slu_pos) + " iso "+str(self.closeness)+"  #Solvate the complex with a cubic solvent box\n")
        
        if self.charge != 0:
            if self.charge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            
            ionnum = str(abs(self.charge))
            ofile.write("addIons2 mol " + ion + ' ' + ionnum + '\n')
        ofile.write('check mol\n')
        
        ofile.write('savepdb mol '+os.path.splitext(self.xyzfile)[0]+'_solvated.pdb\n')
        ofile.write('saveamberparm mol '+os.path.splitext(self.xyzfile)[0]+'_solvated.prmtop' + ' ' + os.path.splitext(self.xyzfile)[0]+'_solvated.inpcrd\n')
        ofile.write('quit\n')
        ofile.close()
    
    def processPackmolPDB(self):
        r"""
        Convert file for custom solvents like CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        data=open(os.path.splitext(self.xyzfile)[0]+'_packmol.pdb','r').readlines()
        output=open(os.path.splitext(self.xyzfile)[0] + '_processed.pdb','w')
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
        solutePDB = os.path.splitext(self.xyzfile)[0] + '.pdb'
        if self.solvent not in  custom_solv_dict.keys():
            print("The type of solvent is not implemented yet")
            sys.exit()
        else:
            solvPrefix = custom_solv_dict[self.solvent]
            solvent_pdb = solvPrefix+'.pdb'
            solvent_pdb_origin = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data/', solvPrefix, solvent_pdb))
            subprocess.call(['cp',solvent_pdb_origin,solvent_pdb])
            output_pdb =  os.path.splitext(self.xyzfile)[0] + "_packmol.pdb"
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
            self.processPackmolPDB()
    
    def writeTleapcmd_add_solvent_custom(self):
        inputpdb = os.path.splitext(self.xyzfile)[0] + '_processed.pdb'
        solvPrefix = custom_solv_dict[self.solvent]
     #   solvent_pdb = solvPrefix+'.pdb'

        solvent_frcmod = solvPrefix+'.frcmod'
        solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvent_frcmod))

        solvent_prep = solvPrefix+'.prep'
        solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                os.path.join('data',solvPrefix,solvent_prep))
        
        ofile = open("leap_add_solventbox.cmd","w")
        ofile.write("source leaprc.protein.ff14SB\n")
        ofile.write("source leaprc.gaff\n")
        ofile.write("source leaprc.water.tip3p\n")
     #   ofile.write(str(amber_solv_dict[str(self.solvent)][0]))
        symbol = 'X1'
        hybrid = 'sp3'
        output = f'addAtomTypes {{\n    {{ "{symbol}" "{self.atomname}" "{hybrid}" }}\n}}'
        ofile.write(output)
        ofile.write('\n')
        ofile.write(self.resname +' = loadmol2 '+self.process_mol2+'\n') 
        ofile.write('loadamberparams '+self.frcmod+'\n')
        ofile.write('loadamberparams frcmod.ionslm_126_opc\n')
        ofile.write("loadamberparams " + solvent_frcmod_path + "\n")
        ofile.write("loadAmberPrep " + solvent_prep_path + "\n")
        ofile.write('mol = loadpdb '+inputpdb+'\n')
        if self.charge != 0:
            if self.charge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            
            ionnum = str(abs(self.charge))
            ofile.write("addIons2 mol " + ion + ' ' + ionnum + '\n')
        
        ofile.write('check mol\n')
        ofile.write('savepdb mol '+os.path.splitext(self.xyzfile)[0]+'_solvated.pdb\n')
        ofile.write('saveamberparm mol '+os.path.splitext(self.xyzfile)[0]+'_solvated.prmtop' + ' ' + os.path.splitext(self.xyzfile)[0]+'_solvated.inpcrd\n')
        ofile.write('quit\n')
        ofile.close()
    
    def writeTleapcmd_add_solvent_user(self):
        ofile = open("leap_add_solventbox.cmd","w")
        ofile.write("source leaprc.protein.ff14SB\n")
        ofile.write("source leaprc.gaff\n")
        ofile.write("source leaprc.water.tip3p\n")
        
        symbol = 'X1'
        hybrid = 'sp3'
        output = f'addAtomTypes {{\n    {{ "{symbol}" "{self.atomname}" "{hybrid}" }}\n}}'
        ofile.write(output)
        ofile.write('\n')

        ofile.write(self.resname + '= loadmol2 '+self.process_mol2+'\n') 
        ofile.write('loadamberparams '+self.frcmod+'\n')
    
        ofile.write("loadoff " + self.solvent_off + "\n")
        ofile.write("loadamberparams " + self.solvent_frcmod + "\n")
        ofile.write("loadamberparams "+ self.frcmod+"\n")
        
        ofile.write("mol = loadmol2 "+os.path.splitext(self.xyzfile)[0] + "_processed.mol2\n")
        
        ofile.write("solvatebox mol " + self.solvent + " " 
                    + str(self.slu_pos) + " iso 0.8  #Solvate the complex with a cubic solvent box\n")
    
        if self.charge != 0:
            if self.charge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            ofile.write("addIons2 mol " + ion + " 0\n")
        ofile.write("check mol\n")
        ofile.write("savepdb mol " + os.path.splitext(self.xyzfile)[0] + "_solvated.pdb\n")
        ofile.write("saveamberparm mol " + os.path.splitext(self.xyzfile)[0] + "_solvated.prmtop " + os.path.splitext(self.xyzfile)[0] + "_solvated.inpcrd\n")
        ofile.write("quit\n")
        ofile.close()
    
    def tleapcmd(self):
        cmd = 'tleap -f leap_add_solventbox.cmd'
        subprocess.call(cmd,shell=True)
    
    def build(self):
        if self.solvent in amber_solv_dict.keys():
            self.get_mass()
            self.gen_LJ_for_ions()
            self.write_mol2()
            self.write_frcmod()
            self.process_mol2()
            self.write_frcmod()
            self.writeTleapcmd_add_solvent()
            self.tleapcmd()
        
        elif self.solvent in custom_solv_dict.keys():
            self.get_mass()
            self.gen_LJ_for_ions()
            self.write_mol2()
            self.write_frcmod()
            self.process_mol2()
            self.write_frcmod()
            self.packSLUSLV()
            self.processPackmolPDB()
            self.writeTleapcmd_add_solvent_custom()
            self.tleapcmd()
        
        elif len(self.solvent_frcmod) > 0 and len(self.solvent_off) > 0:
            self.get_mass()
            self.gen_LJ_for_ions()
            self.write_mol2()
            self.write_frcmod()
            self.process_mol2()
            self.write_frcmod()
            self.writeTleapcmd_add_solvent_user()
            self.tleapcmd()
