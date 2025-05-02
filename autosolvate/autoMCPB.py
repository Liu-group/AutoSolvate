from openbabel import openbabel
import subprocess
import numpy as np
from ase.io import read, write
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown
import io
from contextlib import redirect_stdout, redirect_stderr
import networkx

import os
import re
import sys
import getopt
from glob import glob
from rdkit import Chem

atom_zoo = ['H','C','N','O','CL','BR','I','P','S','F']


metals = [
    "Li", "Be", "Na", "Mg", "Al", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", 
    "Co", "Ni", "Cu", "Zn", "Ga", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", 
    "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", 
    "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", 
    "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "Fr", "Ra", "Ac", 
    "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No", 
    "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", 
    "Lv", "Ts"
]


cheat_metal = {1:'Cu',2:'Fe',3:'Fe',4:'Mn'}

valence_electrons = {
    'H': 1,
    'C': 4,
    'N': 5,
    'O': 6,
    'CL': 7,
    'BR': 7,
    'I': 7,
    'P': 5,
    'S': 6
}

valence_electrons_full = {
    'H': 2,
    'C': 8,
    'N': 8,
    'O': 8,
    'CL': 8,
    'BR': 8,
    'I': 8,
    'P': 8,
    'S': 8
}

def check_ligand(sdf):
    with open(sdf,'r') as f:
        ligand_charge = 0
        check = True
        linker_dic = {}
        type_dic = {}
        connected = []
        data = f.readlines()
       # atomnumbers = int(data[3].split()[0])
        atomnumbers = int(data[3][:3])
        for i in range(atomnumbers):
            linker_dic[i] = 0
        #bondsnumber = int(data[3].split()[1])
        
        bondsnumber = int(data[3][3:6])
        
       # print(bondsnumber)
        connectioninfo = data[4+atomnumbers:bondsnumber+4+atomnumbers]
        print()
        
       # print(connectioninfo)
        atomtypes = data[4:4+atomnumbers]
     #   print(atomtypes)
        for sort,line in enumerate(atomtypes):
            if len(line.split()) > 3:
             #   print(line)
                type_dic[sort] = line.split()[3].upper()
            
        for line in connectioninfo:
         #   print(line)
            atomx = int(line[:3]) - 1
            atomy = int(line[3:6]) - 1
            connected.append(set([atomx,atomy]))
            bondorder = int(line[8])
            linker_dic[atomx] =linker_dic[atomx] + bondorder
            linker_dic[atomy] =linker_dic[atomy] + bondorder
        #    print(atomx,atomy,bondorder)
        
        for atom in type_dic:
            atometype = type_dic[atom]
            valence_electron = valence_electrons[atometype]
            longpair = valence_electrons_full[atometype] - (int(valence_electron)+ linker_dic[atom])
            if longpair < 0:
                check = False
                print('Warning: Explicit valence for atom ',atom+1, atometype, 'is greater than permitted')

                atoms_linked = []
                atoms_linked_type = []
                for linking in connected:
                    if atom in linking:
                        for i in linking:
                            if i != atom:
                                atoms_linked.append(i+1)
                                atoms_linked_type.append(type_dic[i])
                print('atom',atom + 1, atometype, 'is linked with atom', atoms_linked, atoms_linked_type)
                print('Warning: attempting assign',abs(longpair),'positive charge to',atom+1, atometype)
                ligand_charge = ligand_charge + abs(longpair)
            else:
                ligand_charge = ligand_charge - abs(longpair)
    return check, ligand_charge


def ligandbreakdown_BFGS(xyzfile,metalID):
    prefix = os.path.splitext(xyzfile)[0]
    sdf = prefix + '.sdf'
    print(sdf)
    cmd = 'obabel ' + xyzfile + ' -O ' + sdf
    subprocess.call(cmd,shell=True)
    with open(sdf,'r') as f:
        data = f.readlines()
        linker_dic = {}
        type_dic = {}
        connected = []
        atomnumbers = int(data[3][:3])
        bondsnumber = int(data[3][3:6])
        connectioninfo = data[4+atomnumbers:bondsnumber+4+atomnumbers]
        atomtypes = data[4:4+atomnumbers]
        G = networkx.Graph()
        for line in connectioninfo:
            atomx = int(line[:3]) - 1
            atomy = int(line[3:6]) - 1
            if metalID not in [atomx,atomy]:
                G.add_edge(atomx, atomy)
        molecules = list(networkx.connected_components(G))
        molecules_all = []
        for i in molecules:
            molecules_all.append(list(i))
        return molecules_all
            

class AtomInfo():
    def __init__(self,atomnum,atomtype,coordinate):
        self.atomnum = atomnum
        self.atomtype = atomtype
        self.coordinate = coordinate

class AutoMCPB():
    r"""
    Automatically run MCPB.py for metal-organic complexes.
    
    Parameters
    ----------
    filename : str
        The prefix of ligand-metal complex xyz files (not the whole file name) 
    metal_charge_netcharge: int
        Charge of metal as user set
    mode: str, default: "A"
        "A": automatically assign the charge of each ligand bonded with the metal
        "M": manually assign the charge of each ligand bonded with the metal
    spinmult: int, defaut = 1 or 2 based on the number of electrons (MCPB.py).
        Spinmultiplicity of the ligand-metal complex
    Returns
    -------
    None
        To run AutoMCPB, call build function.

    
    """
    def __init__(self, filename, metal_charge, chargefile,mode,spinmult,round,software,amberhome,cutoff, fakecharge): #liglist='', denticity='',ligcons='',atomsinfo='',ligand_charge='',metal_name=''):
        self.metal_charge = metal_charge
        self.filename = filename
        self.xyzfile = filename + '.xyz'
        self.chargefile = chargefile
        self.mode = mode
        self.spinmult = spinmult
        self.round = round
        self.software = software
        self.amberhome = amberhome
        self.cutoff = cutoff
        self.fakecharge = fakecharge

    def coordinates_reader_xyz(self):
        r"""
        read infos from xyz files
        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        all_atom_infos = []
        with open(self.xyzfile) as f:
            data = f.readlines()
            for ID, line in enumerate(data[2:]):
                if len(line.strip()) > 0:
                    atomtype = str.upper(line.split()[0])
                    x = float(line.split()[1])
                    y = float(line.split()[2])
                    z = float(line.split()[3])
                    atom_info = AtomInfo(ID,atomtype,[x,y,z])
                    all_atom_infos.append(atom_info)
        self.atomsinfo = all_atom_infos
    
    def metal_list(self):
        r"""
        find the atom  of metal
        Parameters
        ----------
        None   
 
        Returns
        self.metal_ID
        self.metal_name
        -------
        None
        """
        complex_mol = mol3D()
        complex_mol.readfromxyz(self.xyzfile)
        metal_list = complex_mol.findMetal(transition_metals_only=False)
        self.cheat_metal = None
        if len(metal_list) > 1:
            raise TypeError("Error: autosolvate can't deal with multi-metal systems")
            sys.exit()
        elif len(metal_list) == 0:
            print('Warning: There is no metal in the xyz file or molsimplify can not support this metal.')
            metalall = []
            with open(self.xyzfile) as f:
                data = f.readlines()
                for ID, line in enumerate(data[2:]):
                    if len(line.strip()) > 0:
                        atomtype = str.capitalize(line.split()[0])
                        if atomtype in metals:
                            metalall.append(atomtype)
            if len(metalall) == 0:
                raise TypeError('Error: There is no metal in the xyz file')
                sys.exit()
            elif len(metalall) > 1:
                raise TypeError("Error: autosolvate can't deal with multi-metal systems")
                sys.exit()
            elif len(metalall) == 1:
                print('Find one metal atom',metalall[0])
                self.metal_name = metalall[0].upper()
                self.cheat_metal = metalall[0].upper()
                prefix = os.path.splitext(self.xyzfile)[0]
                print('Use', cheat_metal[int(self.metal_charge)], 'instead of', metalall[0],'to generate ligand files')
                ### generate cheat_metal_xyz ###
                cheat_metalxyz = open(prefix + '_cheat.xyz','w')
                cheat_metalname = cheat_metal[int(self.metal_charge)]
                with open(self.xyzfile,'r') as f:
                    data = f.readlines()
                    cheat_metalxyz.write(data[0])
                    cheat_metalxyz.write(data[1])
                    for line in data[2:]:
                        if len(line.strip()) > 0:
                            atomtype = str.capitalize(line.split()[0]).strip()
                            if atomtype.upper() == self.metal_name:
                                new_string = line.replace(atomtype,cheat_metalname.capitalize())
                            #    print(new_string)
                                cheat_metalxyz.write(new_string)
                            else:
                                cheat_metalxyz.write(line)
                
                cheat_metalxyz.close()
                self.xyzfile = prefix + '_cheat.xyz'
                print('now,the xyzfile is ', self.xyzfile)
                cheat_complex_mol = mol3D()
            #    print(self.xyzfile)
                cheat_complex_mol.readfromxyz(self.xyzfile)
                metal_list = cheat_complex_mol.findMetal(transition_metals_only=False)
                print(metal_list)
                self.metal_ID = metal_list[0]
                self.metal_name =  cheat_metalname.upper()
              #  print(self.metal_name)
                self.coordinates_reader_xyz()
                
        elif len(metal_list) == 1:
            self.cheat_metal = None
            self.metal_ID = metal_list[0]
            metal_type = self.atomsinfo[self.metal_ID].atomtype
            self.metal_name = metal_type
    
    def check_atoms(self):
        r"""
        check the elements whether in atom type of GAFF
        Parameters
        ----------
        None   
        Returns
        -------
        None
        """
        check = 'Yes'
        for atom in self.atomsinfo:
            if atom.atomtype not in atom_zoo + [self.metal_name] :
                print('Error: autosolvate does not support the atom ' + atom.atomtype)
                check = 'No'
                sys.exit()
        
    def ligand_identify(self):
        r"""
        parameters
        ----------
        None

        Returns
        self.liglist: list of list of int, atoms of each ligands [[LG1],[LG2],[LG3]]
        self.ligcons: list of int, atom numbers in complex xyz file, which are bonded with metal
        -------
        None
        """
        complex_mol = mol3D()
        complex_mol.readfromxyz(self.xyzfile)
        ligands_breakdown_infos = ligand_breakdown(complex_mol,transition_metals_only=False)
        liglist = ligands_breakdown_infos[0]
       # denticity = ligands_breakdown_infos[1]
        ligcons =  ligands_breakdown_infos[2]
        if len(liglist) > 0:
            self.liglist = liglist
            self.ligcons = ligcons
        #self.denticity = denticity
        elif len(liglist) == 0:
            print('Warning: No ligand is identified by molsimplify, try to use graph search to break ligand')
            liglist = ligandbreakdown_BFGS(self.xyzfile,int(self.metal_ID))
            self.ligcons = 'NA'
            self.liglist = liglist
            
            
    def modifiy_pdb(self,inputpdb,mol2,outputpdb):
        r"""
        change the format of pdb files to be suitable for GAFF
        parameters
        ----------
        inputpdb: str, .pdb file of ligand, resi ID of this pdb file must be same as the file prefix
        mol2: .mol2 file of ligand
        outputpdb: str, .pdb file name of ligand, LG1.pdb, LG2.pdb ..., 

        Returns
        -------
        None
        """
        inputpdb_atoms = []
        with open(inputpdb,'r') as inf:
            ifile = inf.readlines()
            for line in ifile:
                if line.split()[0] in ['ATOM','HETATM']:
                    if len(line.split()) > 3:
                        inputpdb_atoms.append(line)
     
        with open(mol2,'r') as f:
            data = f.read()
            match = re.search(r"@<TRIPOS>ATOM\n(.*?)\n@<TRIPOS>BOND", data, re.DOTALL)
            if match:
                extracted_data = match.group(1).split('\n')
                with open(outputpdb,'w') as outf:
                    for ID, line in enumerate(extracted_data):
                        new_atomname = line.split()[1]
                        oldline = inputpdb_atoms[ID]
                        newline = oldline[:12] + new_atomname.ljust(4) + oldline[16:]
                        outf.write(newline)
   
    def check_radical(self, ligandxyz): #### this code is for checking the radical linkage
        r"""
        check the radical state of ligand.xyz file to determine the charge that the ligand carry
        parameters
        ----------
        inputpdb: str, ligand.xyz file, like LG1.xyz, LG2.xyz

        Returns
        -------
        charge_of_ligand: int, the charge the ligand carries
        extra_charge: str, 'Y': some extra radical in other atoms which are not bonded with metal
                           'N:: no extra radicals                   
        
        """
        atoms_with_radical = {}
        charge_of_ligand = 0
        
            
        sdf = ligandxyz.split('.xyz')[0]
        ofile = open(sdf + '_' + self.metal_name + '.xyz','w')
        cmd = 'obabel ' + ligandxyz + ' -O ' + sdf + '.sdf'
        with open(sdf + '_sdf.log', 'w') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
    
        heavyatom_of_xyz = []
        checkvalence,attemp_charge = check_ligand(sdf + '.sdf')
        if checkvalence:
            with open(ligandxyz,'r') as f:
                data = f.readlines()
                atomsum_of_ligandxyz = int(data[0])
                atomsum_of_ligandxyz_metal = atomsum_of_ligandxyz + 1
                ofile.write(str(atomsum_of_ligandxyz_metal)+ '\n')
                ofile.write(data[1])
                data_xyz = data[2:]
                for ID,line in enumerate(data_xyz):
                    ofile.write(line)
                    atomtype = line.split()[0]
                    x = float(line.split()[1])
                    y = float(line.split()[2])
                    z = float(line.split()[3])
                    if atomtype != 'H':
                        heavyatom_of_xyz.append(AtomInfo(ID,atomtype,[x,y,z]))
                if len(self.metal_name) == 1:
                    metal_name = self.metal_name
                elif len(self.metal_name) == 2:
                    metal_name = self.metal_name[0] + str.lower(self.metal_name[1])
                else:
                    pass
                 #   print('ERROR: please check the metal_name')
                metal_x = self.atomsinfo[self.metal_ID].coordinate[0]
                metal_y = self.atomsinfo[self.metal_ID].coordinate[1]
                metal_z = self.atomsinfo[self.metal_ID].coordinate[2]
                ofile.write(metal_name+ ' ' + str(metal_x) + ' ' + str(metal_y) + ' ' + str(metal_z))
            ofile.close()

            suppl = Chem.SDMolSupplier(sdf+'.sdf')
            total_radical_number = 0
            mol = [mol for mol in suppl][0]
            
            if mol == None:
                raise TypeError('Error: RDKit can not read the ligand, please check the ligand',ligandxyz)
                sys.exit()
            else:
             #   print('In ',ligandxyz.split('.xyz')[0],':')
                for newID,atom in enumerate(mol.GetAtoms()):
                    radical = atom.GetNumRadicalElectrons()
                    if radical > 0:
                        atoms_with_radical[heavyatom_of_xyz[newID].atomnum] = radical
                        total_radical_number = total_radical_number + radical
                ligand_metal = mol3D()
                ligand_metal.readfromxyz(sdf + '_' + self.metal_name + '.xyz')
                ligand_metal_breakdown_infos = ligand_breakdown(ligand_metal,transition_metals_only=False)

                charge_of_ligand =total_radical_number*(-1)
        else:
            print('Warning: please check the ligand',sdf, 'and use option --mode M to set the charge of each ligand in charge.in')
            charge_of_ligand = attemp_charge

     #   print(sdf,'charge is',charge_of_ligand )
        
        return charge_of_ligand

    def ligand_breakdown(self):
        r"""
        breakdown the ligands from 
        parameters
        ----------
        None

        Returns             
        -------
        self.ligand_charge_dic, dictionary, key: ligand name like: 'LG1', 'LG2', self.ligand_charge_dic['LG1'] = charge, int
        """
        ligand_charge_dic = {}
        extra_charges = []
        for ID, ligand in enumerate(self.liglist):
            length = len(ligand)
            ligandname = 'LG' + str(ID)
            with open(ligandname +'.xyz','w') as f:
                f.write(str(length)+'\n')
                f.write(ligandname+ '\n')
                for sort in ligand:
                    atoms = self.atomsinfo[sort]
                    x = atoms.coordinate[0]
                    y = atoms.coordinate[1]
                    z = atoms.coordinate[2]
                    atomtype = atoms.atomtype
                    if len(atomtype) == 2:
                        atomtype = atoms.atomtype[0].upper() + atoms.atomtype[1].lower()
                    f.write(atomtype + ' ' + str(x) + ' ' + str(y) + ' ' + str(z) + '\n')
            
            charge= self.check_radical(ligandname +'.xyz')
            ligand_charge_dic[ligandname] = charge
            cmd = 'obabel ' + ligandname + '.xyz' + ' -O ' + ligandname + '.smi'
            with open(ligandname + '_obabel_smi.log', 'w') as f:
                subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)

            write(ligandname + '_.pdb', read(ligandname + '.xyz'))   ### ASE will recongnized more elements

            cmd = 'cat ' + ligandname + '_.pdb | grep MOL > ' + ligandname + '__.pdb'
            subprocess.call(cmd,shell=True)

            with open(ligandname + '__.pdb','r') as f:
                data = f.read()
                new_data = data.replace('MOL',ligandname)
                with open(ligandname + '___.pdb','w') as f:
                    f.write(new_data)

        self.ligand_charge_dic = ligand_charge_dic
    
    def write_metal_xyz(self):
        r"""
        generate the mol2 file of metal atom
        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        metalname = self.atomsinfo[self.metal_ID].atomtype
        metal_x = self.atomsinfo[self.metal_ID].coordinate[0]
        metal_y = self.atomsinfo[self.metal_ID].coordinate[1]
        metal_z = self.atomsinfo[self.metal_ID].coordinate[2]
        with open(metalname + '.xyz','w') as f:
            f.write('1\n')
            f.write(str(metalname) + '\n')
            f.write(str(metalname) + ' ' + str(metal_x)+ ' ' +str(metal_y)+' '+str(metal_z))
        
        conv = openbabel.OBConversion()
        conv.SetInAndOutFormats("xyz", "pdb")
        mol = openbabel.OBMol()
        conv.ReadFile(mol, metalname + '.xyz')
        conv.WriteFile(mol, metalname + '_temp.pdb')
        with open(metalname + '_temp.pdb','r') as f:
            data = f.readlines()
        with open(metalname + '.pdb','w') as outf:
            for line in data:
                if 'HETATM' in line:
                    newline = line.replace('UNL',metalname + ' ')
                    new_newline = newline[:12] + str.upper(metalname).ljust(4) + newline[16:76] + str.upper(metalname).rjust(2) + '\n'
                    outf.write(new_newline )
        with open('genmetalmol2.py','w') as f:
            cmd =self.amberhome +'metalpdb2mol2.py -i ' + metalname + '.pdb'+ ' -o ' + metalname + '.mol2' + ' -c ' + str(self.metal_charge)
            subprocess.call(cmd, shell=True,stdout=f, stderr=subprocess.STDOUT)
    
    def generate_mol2_by_auto(self):
        for ID, ligand in enumerate(self.liglist):
         #   print(ID,ligand)
            ligandname = 'LG' + str(ID)
            ligand_charge_dic = self.ligand_charge_dic
            charge = ligand_charge_dic[ligandname]
            if len(ligand) > 1:
                if self.fakecharge == 'N':
                    cmd = self.amberhome+'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c bcc -pf y -nc ' + str(charge) + ' -m 2'
                    with open(ligandname + '_antechamber_generate_mol2.log', 'w') as f:
                        subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                    if os.path.exists(ligandname + '_antechamber_generate_mol2.log'):
                        print('antechamber was processed to generate mol2 file, now checking ' + ligandname + '_antechamber_generate_mol2.log')
                        with open(ligandname + '_antechamber_generate_mol2.log','r') as f:
                            for line in f:
                                if 'Cannot properly run' in line:
                                    if 'sqm -O -i sqm.in -o sqm.out' in line:
                                        print('Warning, bcc can not converge, attemping using fake charges to generate mol2')
                                        atom_count = 0
                                        with open(ligandname +'___.pdb','r') as f:
                                            for line in f:
                                                if line.startswith('ATOM'):
                                                    if len(line) > 10:
                                                        atom_count += 1
                                        fake_charge = charge / atom_count
                                        with open(ligandname +'_fake.chg','w') as f:
                                            for i in range(atom_count):
                                                f.write("%.4f" % fake_charge + ' ')
                                        cmd = self.amberhome +'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c rc -cf ' + ligandname +'_fake.chg'
                                        with open(ligandname + '_antechamber_generate_fake_charge_mol2.log', 'w') as f:
                                            subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)

                                if 'Please check the total charge (-nc flag) and spin multiplicity' in line:
                                    print('warning: the provided charge might not be right for -bcc charge calculation, please check, now attempting using fake charges to generate mol2')
                                    atom_count = 0
                                    with open(ligandname +'___.pdb','r') as f:
                                        for line in f:
                                            if line.startswith('ATOM'):
                                                if len(line) > 10:
                                                    atom_count += 1
                                        fake_charge = charge / atom_count
                                    with open(ligandname +'_fake.chg','w') as f:
                                        for i in range(atom_count):
                                            f.write("%.4f" % fake_charge + ' ')
                                    cmd = self.amberhome +'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c rc -cf ' + ligandname +'_fake.chg'
                                    with open(ligandname + '_antechamber_generate_fake_charge_mol2.log', 'w') as f:
                                        subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                elif self.fakecharge == 'Y':
                    atom_count = 0
                    with open(ligandname +'___.pdb','r') as f:
                        for line in f:
                            if line.startswith('ATOM'):
                                if len(line) > 10:
                                    atom_count += 1
                    fake_charge = charge / atom_count
                    with open(ligandname +'_fake.chg','w') as f:
                        for i in range(atom_count):
                            f.write("%.4f" % fake_charge + ' ')
                    cmd = self.amberhome +'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c rc -cf ' + ligandname +'_fake.chg'
                    with open(ligandname + '_antechamber_generate_fake_charge_mol2.log', 'w') as f:
                        subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                    
                                        
            elif len(ligand) == 1:
                cmd = self.amberhome+'metalpdb2mol2.py -i ' + ligandname + '___.pdb'+ ' -o ' + ligandname + '.mol2' + ' -c ' + str(charge)
                with open(ligandname + '_metalpdb2mol2_generate_mol2.log', 'w') as f:
                    subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                        
    def generate_mol2_by_manual(self):
        r"""
        manually generate the mol2 file if self.mode = 'M'
        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        ligand_length_dic = {}
        for ID, ligand in enumerate(self.liglist):
            ligandname = 'LG' + str(ID)
            length_ligand = len(ligand)
            ligand_length_dic[ligandname] = length_ligand

        ligand_charge_dic = {}
        if self.chargefile == None:
            print('Error: please give a charge input file of charges of all ligands')
        else:
            with open(self.chargefile) as f:
                data = f.readlines()
                if len(data) != len(self.liglist):
                    print('Error: please give the right format of charge input file.')
                for line in data:
                    if len(line.split()) != 2:
                        print('Error: please give the right format of charge input file.')
                    else:
                        ligandname = line.split()[0]
                        ligandcharge = float(line.split()[1])
                        ligand_charge_dic[ligandname] = ligandcharge

        for ligandname in ligand_charge_dic:
            charge = ligand_charge_dic[ligandname]
            length = ligand_length_dic[ligandname]
            if length > 1:
                cmd = self.amberhome+'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c bcc -pf y -nc ' + str(charge) + ' -m 2'
                with open(ligandname + '_antechamber_generate_mol2.log', 'w') as f:
                    subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                
                if os.path.exists(ligandname + '_antechamber_generate_mol2.log'):
                    print('antechamber was processed to generate mol2 file, now checking ' + ligandname + '_antechamber_generate_mol2.log')
                    with open(ligandname + '_antechamber_generate_mol2.log','r') as f:
                        for line in f:
                            if 'Cannot properly run' in line:
                                if 'sqm -O -i sqm.in -o sqm.out' in line:
                                    print('Warning, bcc can not converge, attemping using fake charges to generate mol2')
                                    atom_count = 0
                                    with open(ligandname +'___.pdb','r') as f:
                                        for line in f:
                                            if line.startswith('ATOM'):
                                                if len(line) > 10:
                                                    atom_count += 1
                                    fake_charge = charge / atom_count
                                    with open(ligandname +'_fake.chg','w') as f:
                                        for i in range(atom_count):
                                            f.write("%.4f" % fake_charge + ' ')
                                    cmd = self.amberhome +'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c rc -cf ' + ligandname +'_fake.chg'
                                    with open(ligandname + '_antechamber_generate_fake_charge_mol2.log', 'w') as f:
                                        subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
                            
                            if 'Please check the total charge (-nc flag) and spin multiplicity' in line:
                                print('warning: the provided charge might not be right for -bcc charge calculation, please check, now attempting using fake charges to generate mol2')
                                atom_count = 0
                                with open(ligandname +'___.pdb','r') as f:
                                    for line in f:
                                        if line.startswith('ATOM'):
                                            if len(line) > 10:
                                                atom_count += 1
                                    fake_charge = charge / atom_count
                                with open(ligandname +'_fake.chg','w') as f:
                                    for i in range(atom_count):
                                        f.write("%.4f" % fake_charge + ' ')
                                cmd = self.amberhome +'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c rc -cf ' + ligandname +'_fake.chg'
                                with open(ligandname + '_antechamber_generate_fake_charge_mol2.log', 'w') as f:
                                    subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
            
            elif length == 1:
                cmd = self.amberhome+'metalpdb2mol2.py -i ' + ligandname + '___.pdb'+ ' -o ' + ligandname + '.mol2' + ' -c ' + str(charge)
                with open(ligandname + '_metalpdb2mol2_generate_mol2.log', 'w') as f:
                    subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
        self.ligand_charge_dic = ligand_charge_dic
    
    def rewrite_pdb(self):
        r"""
        generate the frcmod files of ligands
        modify the ligand pdb file        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        for ligandname in self.ligand_charge_dic:
            if os.path.exists(ligandname + '.mol2'):
                cmd = self.amberhome + 'parmchk2 -i ' + ligandname + '.mol2' + ' -o ' + ligandname +'.frcmod' +' -f mol2'
                subprocess.call(cmd,shell=True)
                self.modifiy_pdb(ligandname +'___.pdb',ligandname + '.mol2',ligandname  +'_temp.pdb')
                cmd = 'cat '+ ligandname  + '_temp.pdb  | grep ' + ligandname + ' > ' + ligandname  + '.pdb'  
                subprocess.call(cmd,shell=True)   
                
            else:
                print('Error: ' + ligandname + '.mol2', 'file', 'is not generated by antechamber, please check', ligandname + '_antechamber_generate_mol2.log')
                sys.exit()

    def combine_pdbs(self):
        r"""
        combine the modified pdb files of ligands into a final pdb file
        ----------
        None   
 
        Returns
        -------
        None
        """
        temppdb = self.filename + '_temp.pdb'
        metalname = self.atomsinfo[self.metal_ID].atomtype
        cmd = 'cat ' + metalname + '.pdb '
        for ID, ligandlist in enumerate(self.liglist):
            ligandname = 'LG' + str(ID)
            cmd = cmd + ' ' + ligandname + '.pdb '
        cmd = cmd + '> ' + temppdb
        print('combine ' + cmd)
        subprocess.call(cmd,shell=True)
        cmd =self.amberhome + 'pdb4amber -i ' + temppdb + ' -o ' +  self.filename + '_final.pdb'
        with open(temppdb + '_pdb4amber.log', 'w') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)

    def get_bonded_pairs(self):
        r"""
        get the bonded_pair between 
        ----------
        None   
 
        Returns
        -------
        None
        """
        if self.ligcons != 'NA':
            pdbin = self.filename + '_final.pdb'
            xyzout = self.filename + '_final.xyz'
            write(xyzout, read(pdbin))        
            new_complex_mol = mol3D()
            new_complex_mol.readfromxyz(self.filename + '_final.xyz',)
            metalID = new_complex_mol.findMetal(transition_metals_only=False)
            self.new_complex_mol = new_complex_mol
            pairs = ligand_breakdown(new_complex_mol,transition_metals_only=False)[2]

            add_bonded_pairs = 'add_bonded_pairs '
            for denticitys in pairs:
                for denticity in denticitys:
                    metal_denticity = str(metalID[0]+1) + '-' + str(denticity+1) + ' '
                    add_bonded_pairs = add_bonded_pairs + metal_denticity
            
            self.add_bonded_pairs = add_bonded_pairs
        else:
            self.add_bonded_pairs = ''
        
    def generate_MCPB_input(self):
        metalname = self.atomsinfo[self.metal_ID].atomtype
        if self.ligcons == 'NA':
            pdbin = self.filename + '_final.pdb'
            xyzout = self.filename + '_final.xyz'
            write(xyzout, read(pdbin))        
            new_complex_mol = mol3D()
            new_complex_mol.readfromxyz(self.filename + '_final.xyz',)
            metalID = new_complex_mol.findMetal(transition_metals_only=False)
            self.new_complex_mol = new_complex_mol
        metalID = self.new_complex_mol.findMetal(transition_metals_only=False)
        ion_ids = 'ion_ids ' + str(metalID[0]+1)
        self.ion_ids = ion_ids
        original_pdb = 'original_pdb '+ self.filename + '_final.pdb'
        ion_mol2files = 'ion_mol2files ' + metalname +'.mol2'
        naa_mol2files = 'naa_mol2files '
        frcmod_files = 'frcmod_files '
        for ID, ligandlist in enumerate(self.liglist):
            ligandname = 'LG' + str(ID)
            lgmol2 = ligandname + '.mol2'
            lgfrcmod = ligandname + '.frcmod'
            naa_mol2files = naa_mol2files + ' ' + lgmol2
            frcmod_files = frcmod_files + ' ' + lgfrcmod
        self.naa_mol2files = naa_mol2files
        self.frcmod_files = frcmod_files

        with open(self.filename + '_MCPB.in','w') as f:
            f.write('software_version ' + self.software +'\n')
            f.write(original_pdb + '\n')
            f.write('group_name ' + self.filename + '\n')
            f.write('cut_off ' + str(self.cutoff) + '\n')
            f.write(ion_ids + '\n')
            if self.ligcons != 'NA':
                f.write(self.add_bonded_pairs + '\n')
            f.write(ion_mol2files+ '\n')
            f.write(naa_mol2files+ '\n')
            f.write(frcmod_files + '\n')
            f.write('force_field  ff14SB\n')
            f.write('gaff 1\n')
            if self.spinmult == None:
                pass
               # print('the spin is default by the number of electrons')
            else:
              #  print('the spin is set as '+ str(self.spinmult))
                f.write('lgmodel_spin ' + str(self.spinmult)+'\n')
                f.write('smmodel_spin ' + str(self.spinmult)+'\n')
            
    
    def find_atomnumber(self, complex_xyz, ligand_metal_xyz):
        epsilon = 1e-3
        atom_number_dic = {}
        ligand_metal_atom_info = []
        complex_atom_info = []
        with open(ligand_metal_xyz,'r') as f:
            data = f.readlines()
            for ID,line in enumerate(data[2:]):
                atomtype = line.split()[0]
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                ligand_metal_atom_info.append(AtomInfo(atomtype=atomtype,coordinate=[x,y,z],atomnum=ID))
        
        with open(complex_xyz,'r') as f:
            data = f.readlines()
            for ID, line in enumerate(data[2:]):
                atomtype = line.split()[0]
                x = float(line.split()[1])
                y = float(line.split()[2])
                z = float(line.split()[3])
                complex_atom_info.append(AtomInfo(atomtype=atomtype,coordinate=[x,y,z],atomnum=ID))

        for atom in ligand_metal_atom_info:
            atomnumber = atom.atomnum
            coordinate = np.array(atom.coordinate)
            atomtype  = atom.atomtype
            for item in complex_atom_info:
                if str.upper(atomtype) == str.upper(item.atomtype):
                    are_close = np.allclose(coordinate, np.array(item.coordinate),atol=epsilon)
                    if are_close:
                        atom_number_dic[atomnumber] = item.atomnum
        return(atom_number_dic)
    
    def get_connection_info(self): ### to write the charge and connection bond type, and pair of connection
        with open(self.filename + '.info','w') as f:
            f.write('@charges\n')
            f.write(str.upper(self.metal_name) + ' ' + str(self.metal_charge) + '\n')
            for ligand in self.ligand_charge_dic:
                f.write(ligand + ' ' + str(self.ligand_charge_dic[ligand]) + '\n')
            
            f.write('@connection_infos\n')
        #   for ligand in self.ligand_charge_dic:
         #       connections = self.coorination_or_covalent(ligand, self.filename + '_final.xyz', ligand + '_' + str.upper(self.metal_name) + '.xyz')
         #       for connections_info in connections:
         #           f.write(ligand + ' ' + connections_info + '\n')
    
    def get_totalcharge(self):
        totalcharge = 0
        for ligand in self.ligand_charge_dic:
            totalcharge = totalcharge + self.ligand_charge_dic[ligand]
        totalcharge = totalcharge + self.metal_charge
        return(totalcharge)
    
    def check_MCPBinp(self):
            mcpbinp = self.filename + '_MCPB.in'
            if os.path.exists(mcpbinp):
                with open(mcpbinp,'r') as f:
                    for line in f:
                        if 'naa_mol2files' in line:
                            mol2file = line.split()[1:]
                            for mol2 in mol2file:
                                if os.path.exists(mol2):
                                    print(mol2, 'is generated')
                                else:
                                    raise TypeError('Erorr: can not find',mol2,'for MCPB input file')
                                    sys.exit()
                        elif "ion_mol2files" in line:
                            mol2file = line.split()[1:]
                            for mol2 in mol2file:
                                if os.path.exists(mol2):
                                    print(mol2, 'is generated')
                                else:
                                    raise TypeError('Erorr: can not find',mol2,'for MCPB input file')
                                    sys.exit()                           
                        elif 'frcmod_files' in line:
                            frcmodfiles = line.split()[1:]
                            for frcmod in frcmodfiles:
                                if os.path.exists(frcmod):
                                    print(frcmod, 'is generated')
                                else:
                                    raise TypeError('Erorr: can not find',frcmod,'for MCPB input file')
                                    sys.exit()                           
                        elif "original_pdb" in line:
                            pdbfile = line.split()[1:]
                            for pdb in pdbfile:
                                if os.path.exists(pdb):
                                    print(pdb,'is generated')
                                else:
                                    raise TypeError('Erorr: can not find',pdb,'for MCPB input file')
                                    sys.exit()     
                        elif "add_bonded_pairs" in line:
                            print(line.strip())
            else:
                raise TypeError('Erorr: can not find',mcpbinp,' for MCPB.py -s 1!')
                sys.exit()
            mcpbinfo =       self.filename + '.info'
            if os.path.exists(mcpbinfo):
                print('charge assigned for each ligand:')
                with open(mcpbinfo,'r') as f:
                    data = f.readlines()
                    for line in data[1:-1]:
                        print(line.strip())
            else:
                raise TypeError('Erorr, can not find',mcpbinfo)
    
    def change_thecheat_atom(self):
        if self.cheat_metal != None:
            prefix = os.path.splitext(self.xyzfile)[0].split('_cheat')[0]
            ##generate real mol2 file##
            newpdb = open(self.cheat_metal.upper()+'.pdb','w')
            with open(self.metal_name.upper() + '.pdb','r') as f:
                for line in f:
                    if 'ATOM' or 'HETATM' in line.split()[0]:
                        newline = line.replace(self.metal_name.upper(),self.cheat_metal.upper())
                        newpdb.write(newline)
            newpdb.close()
            cmd = 'rm ' + self.metal_name.upper() + '.pdb'
            subprocess.call(cmd,shell=True)
            
            cmd = 'rm ' + self.metal_name.upper() + '_temp.pdb'
            subprocess.call(cmd,shell=True)           
        
            cmd =self.amberhome +'metalpdb2mol2.py -i ' + self.cheat_metal.upper()+'.pdb' + ' -o ' + self.cheat_metal.upper() + '.mol2' + ' -c ' + str(self.metal_charge)
            subprocess.call(cmd, shell=True)
            
            cmd = 'rm ' + self.metal_name.upper() + '.mol2'
            subprocess.call(cmd,shell=True)
            
            ### generate real final.pdb ###
            finalpdbname = prefix + '_final.pdb'
            cmd = 'cp ' + finalpdbname + ' ' + prefix + '_final_old.pdb'
            subprocess.call(cmd,shell=True)
            
            finalpdb = open(finalpdbname,'w')
            with open(prefix + '_final_old.pdb','r') as f:
                data = f.readlines()
                for line in data:
                    if 'ATOM' or 'HETATM' in line.split()[0]:
                        if self.metal_name.upper() in line[12:16]:
                            newline = line[:6] + line[6:].replace(self.metal_name.upper(),self.cheat_metal.upper())
                            finalpdb.write(newline)
                        else:
                            finalpdb.write(line)
            finalpdb.close()
            
            pdbin = prefix + '_final.pdb'
            xyzout = prefix + '_final.xyz'
            write(xyzout, read(pdbin))   
            
            ### generate real info file ###
            
            infofile = prefix + '.info'
            cmd = 'cp ' + infofile  + ' ' + prefix + '_old.info'
            subprocess.call(cmd,shell=True)
            
            finalinfo = open(infofile,'w')
            with open(prefix + '_old.info','r') as f:
                data = f.readlines()
                for line in data:
                    if self.metal_name.upper() in line:
                        newline = line.replace(self.metal_name.upper(),self.cheat_metal.upper())
                        finalinfo.write(newline)
                    else:
                        finalinfo.write(line)
            finalinfo.close()
            
            ### generate real MCPB input ###
            MCPBin = prefix + '_MCPB.in'
            cmd = 'cp ' + MCPBin + ' ' + prefix + '_MCPB_old.in'
            subprocess.call(cmd,shell=True)
            
            finalMCPBinput = open(MCPBin,'w')
            with open(prefix + '_MCPB_old.in','r') as f:
                data = f.readlines()
                for line in data:
                    if 'ion_mol2files' in line:
                        newline = 'ion_mol2files ' + self.cheat_metal.upper() + '.mol2\n'
                        finalMCPBinput.write(newline)
                    else:
                        finalMCPBinput.write(line)
            finalMCPBinput.close()
            
            if self.software == 'orca':
                MCPBin = prefix + '_MCPB_orca.in'
                cmd = 'cp ' + MCPBin + ' ' + prefix + '_MCPB_orca_old.in'
                subprocess.call(cmd,shell=True)
                
                finalMCPBinput = open(MCPBin,'w')
                with open(prefix + '_MCPB_old.in','r') as f:
                    data = f.readlines()
                    for line in data:
                        if 'ion_mol2files' in line:
                            newline = 'ion_mol2files ' + self.cheat_metal.upper() + '.mol2\n'
                            finalMCPBinput.write(newline)
                        else:
                            finalMCPBinput.write(line)
                finalMCPBinput.close()
            
            
    
    def run_MCPB_1(self):
        if self.round in ['1']:
            cmd =self.amberhome +'MCPB.py -i ' + self.filename + '_MCPB.in -s ' + self.round
            if self.software in ['orca','ORCA','Orca']:
                with open(self.filename + '_MCPB.in','r') as source_file:
                    data = source_file.readlines()
                with open(self.filename + '_MCPB_orca.in','w') as target_file:
                    target_file.write('software_version gms\n')
                    for line in data:
                        if 'software_version orca' not in line:
                            target_file.write(line)
                cmd = self.amberhome +'MCPB.py -i ' + self.filename + '_MCPB_orca.in -s ' + self.round
            with open('MCPB_1.log', 'w') as output_file:
                subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
            
        #    totalcharge = self.get_totalcharge()
    
    def run_MCPB_2(self):    
        if  self.round in ['2','2b','2e','2s','2z']:
          #  print('******** start secondary round of MCPB.py ********\n' )
            print('checking the file generated by QM calculation:\n')
            if self.software in ['g09','gau','g03']:
                with open('MCPB_2.log','w') as f:
                    cmd = self.amberhome +'MCPB.py -i ' + self.filename + '_MCPB.in -s ' + self.round
                    subprocess.call(cmd,shell=True,stdout=f,stderr=subprocess.STDOUT)
            elif self.software in ['gms','GAMESS','GMS']:
                for resultfile in [self.filename + '_small_fc.log',self.filename + '_small_opt.log',self.filename + '_large_mk.log']:
                    if  os.path.exists(resultfile):
                        print(resultfile  + ' is generated by the first round of MCPB.py \n' )
                    else:
                        print('Error', resultfile  + ' is not generated by the first round of MCPB.py \n' )
                        sys.exit()
                cmd = self.amberhome + 'MCPB.py -i ' + self.filename + '_MCPB.in -s ' + self.round
                subprocess.call(cmd,shell=True)
            elif self.software in ['orca','ORCA','Orca']:
                for resultfile in [self.filename + '_small_fc.orcaout',self.filename + '_small_opt.orca'+'_trj.xyz',
                                   self.filename + '_small_opt.orcaout',self.filename + '_large_mk.orcaout']:
                    if os.path.exists(resultfile) :
                        print(resultfile  + ' is generated by the first round of MCPB.py \n' )
                    else:
                        print('Error', resultfile  + ' is not generated by the first round of MCPB.py \n' )
                        sys.exit()
                cmd = 'MCPB_orca ' + self.filename + '_MCPB.in ' + self.round
                with open('MCPB_2.log', 'w') as output_file:
                    subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)


    
    def run_MCPB_3(self): 
        if self.round in ['3','3b','3c','3d']:
          #  print('******** start third round of MCPB.py ********')
            for resultfile in [self.filename + '_mcpbpy.frcmod',self.filename + '_mcpbpy_pre.frcmod']:
                if os.path.exists(resultfile):
                    pass
                 #   print(resultfile  + ' is generated by the secondary round of MCPB.py' )
                else:
                    print('Error', resultfile  + ' is not generated by the secondary round of MCPB.py')
            cmd =self.amberhome +'MCPB.py -i ' + self.filename + '_MCPB.in -s ' + self.round
            if self.software in ['orca','ORCA','Orca']:
                cmd = 'MCPB_orca ' + self.filename + '_MCPB.in ' + self.round
            with open('MCPB_3.log', 'w') as output_file:
                subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
    
    def check_tleapin(self, inp, outfile):
        ifile = open(inp, 'r')
        data_tleapin = ifile.readlines()
        ifile.close()
        out = open(outfile, 'w')
        AtomTypes = []
        bondmol = []  
        
        for line in data_tleapin:
            if '{' in line and '}' in line:
                if line not in AtomTypes:
                    out.write(line)
                    AtomTypes.append(line)
            elif 'bond mol.' in line:
                if line not in bondmol:  
                    out.write(line)
                    bondmol.append(line)  
            else:
                out.write(line)
        out.close()


    def run_MCPB_4(self): 
        if self.round in ['4','4b','4n1','4n2']:
         #   print('******** start fourth round of MCPB.py ********\n' )
            for resultfile in [self.filename + '_mcpbpy.frcmod',self.filename + '_mcpbpy_pre.frcmod']:
                if  os.path.exists(resultfile):
                    pass
                    print(resultfile  + ' is generated by the secondary round of MCPB.py \n' )
                else:
                    print('Error', resultfile  + ' is not generated by the secondary round of MCPB.py \n' )
            cmd = self.amberhome +'MCPB.py -i ' + self.filename + '_MCPB.in -s ' + self.round
            if self.software in ['orca','ORCA','Orca']:
                cmd = 'MCPB_orca ' + self.filename + '_MCPB.in ' + self.round
            with open('MCPB_4.log', 'w') as output_file:
                subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
            inp = self.filename + '_tleap.in'
            outfile = self.filename + '_tleap_check.in'
            self.check_tleapin(inp,outfile)
            cmd =self.amberhome +'tleap -f ' + outfile + ' > ' + outfile.split('.in')[0] + '.out'
            with open('tleap_MCPB.log','w') as output_file: 
                subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)

    def checkingFF(self):
        print('******** start to check the force field ******** ')
        missingbond = open('missingbonds.txt','w')
        topology_file = self.filename + '_dry.prmtop'
        cpptraj_command = 'bondinfo\n'
        with open('bondinfo.in','w') as f:
            f.write(cpptraj_command)
        with open('bondinfo_output.txt', 'w') as output_file:
            cmd = self.amberhome + 'cpptraj -p ' + topology_file + ' -i bondinfo.in'
            subprocess.call(cmd,shell=True,stdout=output_file, stderr=subprocess.STDOUT)
       
        with open(self.filename + '_MCPB.in','r') as file:
            for line in file:
                if 'ion_ids' in line:
                    ion_ids = line.split()[1]
                elif 'add_bonded_pairs' in line:
                    add_bonded_pairs = line.split()[1:]
                  #  print(add_bonded_pairs)
                elif 'ion_mol2files' in line:
                    metal_name = line.split()[1].split('.')[0]
        with open('bondinfo_output.txt','r') as f:
            output = f.read()
        for sort, line in enumerate(output.split('\n')):
            if '#Bnd'  in line:
                begin = sort + 1
            elif 'TIME' in line:
                end = sort
            elif 'Error' in line:
                print('try to resort atom number')
                with open('fixatord.in','w') as cpptrajin:
                    cpptrajin.write('fixatomorder outprefix fixed\n')
                    cpptrajin.write('trajout restart' + 'fixed_' + self.filename + '_dry.inpcrd' )
                    cpptrajin.write('trajout restart' + 'fixed_' + self.filename + '_dry.prmtop' )
                    cpptrajin.write('run\n')
                    cpptrajin.write('quit\n')
        
                cmd = self.amberhome + 'cpptraj -p ' + self.filename + '_dry.prmtop -c ' + self.filename +  '_solv.inpcrd -i fixatord.in > fixatord.out'
                subprocess.call(cmd,shell=True)

        bonded_pairs_in_prmtop = []
        for line in output.split('\n')[begin:end]:
            A1 = line.split()[5]
            A2 = line.split()[6]
            pairlist = [A1,A2]
            if ion_ids in [A1,A2]:
                pairlist.remove(ion_ids)
                bond_pair = ion_ids + '-' + pairlist[0]
                bonded_pairs_in_prmtop.append(bond_pair)
        
        missingbondall = []
        
        for pair in add_bonded_pairs:
            if pair in bonded_pairs_in_prmtop:
                print('metal bond ' + pair + ' is parameterized in prmtop file')
            else:
                print('Warning: metal bond ' + pair + ' is not parameterized in prmtop file')
                missingbondall.append(pair)
        if len(missingbondall) == 0:
            missingbond.write('None')
            missingbond.close()
        else:
            for pair in missingbondall:
                missingbond.write(pair+'\n')
        
        for pair in bonded_pairs_in_prmtop:
            if pair not in add_bonded_pairs:
                print(pair)
                print(add_bonded_pairs)
                print('extra metal bond ' + pair + ' is not assigned in MCPB.in but parameterized in prmtop file')
    
        print('******** start to use ParmEd to check Force Field ********')
        with open('mcpbpy_parmed.in','w') as parmedfile:
            metalID  = metal_name + ion_ids
            parmedfile.write('printBonds :'+ metalID + '\n')
            parmedfile.write('printAngles :'+ metalID + '\n')
            parmedfile.write('printDihedrals :' + metalID + '\n')
            parmedfile.write('printDetails :' + metalID + '\n')
        
        cmd = 'parmed -i mcpbpy_parmed.in -p ' + self.filename + '_dry.prmtop > parmed.out'
        subprocess.call(cmd,shell=True)

        with open('parmed.out','r') as parmedout:
            data = parmedout.readlines()
        ### check bonds###
        bond_K_b = 'NAN'
        bond_K_e = 'NAN'
        angle_b = 'NAN'
        angle_e = 'NAN'
        dihedral_b = 'NAN'
        dihedral_b = 'NAN'
        for sort,line in enumerate(data):
            if 'Atom 1              Atom 2       R eq   Frc Cnst' in line:
                bond_K_b = sort + 1
            elif 'Atom 1               Atom 2               Atom 3   Frc Cnst   Theta eq' in line:
                bond_K_e = sort - 1
                angle_b = sort + 1
            elif 'Phase  EEL Scale  VDW Scale' in line:
                angle_e = sort - 1
                dihedral_b = sort + 1                   
            elif "Charge GB Radius GB Screen" in line:
                FF_of_metal = data[sort + 1]
                dihedral_e = sort - 4
                LJ_Radius = FF_of_metal.split()[6]
                if float(LJ_Radius) > 1:
                    print('LJ_Radius of metal bonds' + ' > 1 Angstrom')
                elif float(LJ_Radius) < 1:
                    print('Warning: LJ__Radius of metal bonds < 1 Angstrom')
                
                LJ_Depth = FF_of_metal.split()[7]
                Mass = FF_of_metal.split()[8]
                Charge = FF_of_metal.split()[9]
                GB_Radius = FF_of_metal.split()[10]
                GB_Screen = FF_of_metal.split()[11]
            
        bond_fcs = data[bond_K_b:bond_K_e]
        bond_fcs_all = []
        bond_R_all = []

        for line in bond_fcs:
            bond_frc= float(line.split()[-1])
            bond_R = float(line.split()[-2])
            bond_fcs_all.append(bond_frc)
            bond_R_all.append(bond_R)
            bond_infos = line.split()[:-2]
            bond_info = ''
            for line in bond_infos:
                bond_info =bond_info + ' ' + line 
            if bond_frc < 200:
                pass
            elif bond_frc >= 200:
                print('Warning bond Frc Cnst:',bond_info, '>=','200')
            if bond_R < 2.8:
                pass
            elif bond_R >= 2.8:
                print('Warning bond R eq:',bond_info, '<','200')
        
        bond_fcs_all_avg = sum(bond_fcs_all)/len(bond_fcs_all)
        bond_R_all_avg = sum(bond_R_all)/len(bond_R_all)
        if bond_fcs_all_avg >= 200:
            print('Warning:', 'the avg of bond Frc Cnst > 200' )
        else:
            print('the avg of bond Frc Cnst < 200 ')
        if bond_R_all_avg >= 2.8:
            print('Warning:', 'the avg of bond R eq > 2.8' )
        else:
            print('the avg of bond R e < 200 ')
        
        if angle_b != 'NAN':
            if angle_e != 'NAN':
                angle_fcs = data[angle_b:angle_e]
                if len(angle_fcs) >= 1:
                    angle_fcs_all = []
                    theta_all = []
                    for line in angle_fcs:
                        angle_fc = float(line.split()[-2])
                        theta = float(line.split()[-1])
                        angle_fcs_all.append(angle_fc)
                        theta_all.append(theta)
                        new_line = line.replace(line.split()[-1]+'\n','')
                        new_new_line = new_line.replace(line.split()[-2],'')
                        if theta <= 100 :
                            pass
                           # print('Warning: Theta eq of ', new_new_line, 'is < 100 degree')
                        if angle_fc >= 100:
                            pass
                         #   print('Warning: Frc Cnst of ', new_new_line, 'is > 100')


                        angle_fc_avg = sum(angle_fcs_all)/len(angle_fcs_all)
                        theta_fc_avg = sum(theta_all)/len(theta_all)
                    if angle_fc_avg > 100:
                        print('Warning: the avg of angle Frc Cnst > 100 ' )
                    else:
                        print('the avg of angle Frc Cnst < 100 ' )
                    
                    if theta_fc_avg >= 100:
                        print('Warning: the avg of  Theta eq >= 100 degree' )
                    else:
                        print('the avg of Theta eq < 100 degree' )

        if dihedral_b != 'NAN':
            if dihedral_e != 'NAN':
                dihedrals = data[dihedral_b:dihedral_e]
                dihedrals_all = []
                if len(dihedrals) >= 1:
                    for line in dihedrals:
                        if len(line.strip()) > 5:
                            height = float(line.split()[-5])
                         #   print(height,line)
                            dihedrals_all.append(height)
                            if height > 0:
                                print('Warning: ' + line[:87], 'height','> 0'  )
                    height_avg = sum(dihedrals_all) / len(dihedrals_all)
                    if height_avg > 0:
                        print('Warning: the avg height of metal dihedral > 0')

    def build(self):
        if self.round == '1':
            self.coordinates_reader_xyz()
            self.metal_list()
            self.check_atoms()
            self.ligand_identify()
            self.ligand_breakdown()
            self.write_metal_xyz()
            if self.mode == 'A':
                self.generate_mol2_by_auto()
            elif self.mode == 'M':
                self.generate_mol2_by_manual()
            self.rewrite_pdb()
            self.combine_pdbs()
            self.get_bonded_pairs()
            self.generate_MCPB_input()
            self.get_connection_info()
            self.check_MCPBinp()
            self.change_thecheat_atom()
            self.run_MCPB_1()

        elif self.round in ['2','2b','2e','2s','2z']:
            self.run_MCPB_2()
        elif self.round in ['3','3a','3b','3c','3d']:
            self.run_MCPB_3()
        elif self.round in ['4','4b','4n1','4n2']:
            self.run_MCPB_4()
        elif self.round in ['5','check','CHECK']:
            self.coordinates_reader_xyz()
            self.metal_list()
            self.check_atoms()
            self.ligand_identify()
            if self.ligcons != 'NA':
                self.checkingFF()

def startautoMCPB(argumentList):
    options = "hn:c:u:m:f:s:x:A:e"
    long_options = ["help",'filename=','metal_charge=','spin=','mode=','chargefile=','round=','software=','amberhome=','cutoff=']
    arguments, values = getopt.getopt(argumentList,options,long_options)
    filename = None
    metal_charge = None
    mode='A'
    spinmult=None
    chargefile = None
    software = 'gms'
    round = '1'
    amberhome = '$AMBERHOME/bin/'
    cutoff = '2.8'
   # print(arguments)
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "--help"):
            message = '''Usage: autoMCPB.py [OPTIONS]
            Options:
                -h, --help              Display this help message.
                -n, --filename          Specify the prefix of metal complex.xyz file.
                -c, --metal_charge      Set the metal charge.
                -u, --spin              Set the spin defaut: 1 or 2.
                -m, --mode              Set the mode (A/M).
                -f, --chargefile        Specify the charge file if --mode is M.
                                        chargefile example:
                                        LG1 -1 # ligand name charge 
                -s, --round             round of MCPB.py
                -x, --software            g09,g16 or gms
                -A, --amberhome         path of AmberTools bin
                -e  --cutoff            cutoff of MCPB 
                '''
            print(message)
            sys.exit()
        elif currentArgument in ("-n","--filename"):
       #     print ("metal_complex filename", currentValue)
            filename = str(currentValue)
        elif currentArgument in ("-c",'--metal_charge'):
            if metal_charge == 'default':
                metal_charge == 2
            metal_charge = int(currentValue)
        elif currentArgument in ('-u','--spin'):
            spinmult=int(currentValue)
        elif currentArgument in ('-m','--mode'):
            mode = str(currentValue)
        elif currentArgument in ('-f','--chargefile'):
            chargefile = str(currentValue)
        #    print('read the charge from', chargefile)
        elif currentArgument in ('-s','--round'):
            round = str(currentValue)
        elif currentArgument in ('-x','--software'):
            software = str(currentValue)
        elif currentArgument in ('-A','--amberhome'):
            amberhome = str(currentValue)
        elif currentArgument in ('-e','--cutoff'):
            cutoff = str(currentValue)



    builder = AutoMCPB(filename=filename,metal_charge=metal_charge, spinmult=spinmult,cutoff=cutoff,
                       mode=mode,chargefile=chargefile,round=round,software=software,amberhome=amberhome)

    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startautoMCPB(argumentList)
