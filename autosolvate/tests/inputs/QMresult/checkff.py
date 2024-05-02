import subprocess
import numpy as np
import parmed

def are_points_close(point1, point2):
    answer = False
    diff = [abs(a - b) for a, b in zip(point1, point2)]
    diff =sum(diff)
    if diff < 0.001:
        answer = True
    return answer

def collectdata(sections,):
    alldata = []
    for line in sections:
        alldata.extend(line.split())
    return alldata

def get_atoms_bond(sections,atomtypedic):
    allbonds = {}
    for group in sections:
        atomA  = int(int(group[0])/3) + 1
        atomB = int(int(group[1])/3) + 1
        atomAtypes = atomtypedic[atomA]
        atomBtypes = atomtypedic[atomB]
        allbonds[(atomA,atomB)] =set([atomAtypes,atomBtypes])
    return allbonds

def get_atoms_angle(sections,atomtypedic):
    allbonds = []
    for group in sections:
        atomA  = int(int(group[0])/3) + 1
        atomB = int(int(group[1])/3) + 1
        atomC = int(int(group[2])/3) + 1
        atomAtypes = atomtypedic[atomA]
        atomBtypes = atomtypedic[atomB]
        atomCtypes = atomtypedic[atomC]
        allbonds.append([atomAtypes,atomBtypes,atomCtypes])
    return allbonds

class AtomInfo():
    def __init__(self,atomnum,atomname,atomtype,resname,resseq,x,y,z):
        self.atomnum = atomnum
        self.atomname = atomname
        self.atomtype = atomtype
        self.resname = resname
        self.resseq = resseq
        self.x = x
        self.y = y
        self.z = z

class checkff:
    def __init__(self,  filename):
        self.filename = filename
        self.pdb = self.filename + '_mcpbpy.pdb'
        
    def readpdb(self):
        pdbinfos = []
        lignameinfos = []
        with open(self.pdb,'r') as f:
            for line in f:
                if line.split()[0] in ['ATOM','HETATM']:
                    atomnum = int(line[6:11])
                    atomname = str(line[12:16].strip())
                    resname = str(line[17:20].strip())
                    resseq = str(line[22:26].strip())
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    if resname not in lignameinfos:
                        lignameinfos.append(resname)
                    i  = AtomInfo(atomname=atomname,atomtype=None,
                                  resname=resname,resseq=resseq,atomnum=atomnum,
                                  x=x,y=y,z=z)
                    pdbinfos.append(i)

        return pdbinfos, lignameinfos
    
    def readmol2(self,lignameinfos,pdbinfos):
        start = False
        newpdbinfos = []
        atomtypedic = {}
        for mol2 in lignameinfos:
            with open(mol2+'.mol2','r') as f:
                b = '@<TRIPOS>ATOM'
                e = '@<TRIPOS>BOND'
                for line in f:
                    if b in line:
                        start = True
                        continue
                    if start:
                        if e in line:
                            break
                        if len(line.split()) > 8:
                            atomname = line.split()[1]
                            atomtype = line.split()[5]
                            resname = line.split()[7]
                            x = float(line.split()[2])
                            y = float(line.split()[3])
                            z = float(line.split()[4])
                            for atom in pdbinfos:
                                if atom.atomname == atomname:
                                    if atom.resname == resname:
                                        if are_points_close([x,y,z],[atom.x,atom.y,atom.z]):
                                            atom.atomtype = atomtype
                                            newpdbinfos.append(atom)
                                            atomtypedic[atom.atomnum] = atomtype
        return newpdbinfos,atomtypedic
    
    def readfrcmod(self):
        bondinfos = {}
        angleinfos = {}
        dihedralinfos = {}
        with open(self.filename + '_mcpbpy.frcmod','r') as f:
            data = f.readlines()
            bond = 0
            angl = 0
            dihe = 0
            impr = 0
            for sort,line in enumerate(data):
                if 'BOND' in line:
                    bond = sort + 1
                elif 'ANGL' in line:
                    angl = sort + 1
                elif 'DIHE' in line:
                    dihe = sort + 1
                elif 'IMPR' in line:
                    impr = sort + 1

            bonds = data[bond:angl-1]
            angles = data[angl:dihe-1]
            dihes = data[dihe:impr-1]
            
            for line in bonds:
                if len(line) > 5:
                    pair = line.split()[0]
                    fc = float(line.split()[1])
                    eq = float(line.split()[2])
                    atomA = pair.split('-')[0]
                    atomB = pair.split('-')[1]
                    pair_set = tuple([atomA,atomB])
                    bondinfos[pair_set] = [fc,eq]

            for line in angles:
                if len(line) > 5:
                    if 'nan' not in line:
                        pair = line.split()[0]
                        fc = float(line.split()[1])
                        eq = float(line.split()[2])
                        atomA = pair.split('-')[0]
                        atomB = pair.split('-')[1]
                        atomC = pair.split('-')[2]
                        pair_set = (atomA,atomB,atomC)
                        angleinfos[pair_set] = [fc,eq]
            
            for line in dihes:
                if len(line) > 5:
                    if 'MCPB.py' in line:
                        pair = line.split()[0]
                        fc = float(line.split()[2])
                        eq = float(line.split()[3])
                        atomA = pair.split('-')[0]
                        atomB = pair.split('-')[1]
                        atomC = pair.split('-')[2]
                        atomD = pair.split('-')[3]

                        angleinfos[(atomA,atomB,atomC,atomD)] = [fc,eq]
            
        return bondinfos, angleinfos, dihedralinfos

    
    def gen_new_prmtop(self,oldprmtop,newprmtop,missingbond,atomtypedic,bondinfos,angleinfos):
        with open(missingbond,'r') as f:
            data = f.read()
            if 'None' in data:
                print('There is no missing bonds of metal, just copy the oldprmtop to newprmtop')
                cmd = 'cp ' + oldprmtop + ' ' + newprmtop
                subprocess.call(cmd,shell=True)
            else:
                with open(oldprmtop,'r') as f:
                    data = f.readlines()
                    bonds_with_h = 0
                    bonds_no_h = 0
                    angle_with_h = 0
                    angle_no_h = 0
                    dih = 0
                    for sort, line in enumerate(data):
                        if 'BONDS_INC_HYDROGEN' in line:
                            bonds_with_h = sort + 2
                        elif 'BONDS_WITHOUT_HYDROGEN' in line:
                            bonds_no_h = sort + 2
                        elif 'ANGLES_INC_HYDROGEN' in line:
                            angle_with_h = sort + 2
                        elif 'ANGLES_WITHOUT_HYDROGEN' in line:
                            angle_no_h = sort + 2
                        elif 'DIHEDRALS_INC_HYDROGEN' in line:
                            dih = sort 
                    
                    bonds_with_h_data = collectdata(data[bonds_with_h:bonds_no_h-2])
                    bonds_without_h_data = collectdata(data[bonds_no_h:angle_with_h-2])
                    angle_with_h_data = collectdata(data[angle_with_h:angle_no_h-2])
                    angle_without_h_data = collectdata(data[angle_no_h:dih])
                    print(angle_without_h_data)

                    groups_bonds_with_h = [bonds_with_h_data[i:i+3] for i in range(0, len(bonds_with_h_data), 3)]
                    groups_bonds_without_h = [bonds_without_h_data[i:i+3] for i in range(0, len(bonds_without_h_data), 3)]
                    groups_angles_with_h =  [angle_with_h_data[i:i+4] for i in range(0, len(angle_with_h_data), 4)]
                    groups_angles_without_h = [angle_without_h_data[i:i+4] for i in range(0, len(angle_without_h_data), 4)]

                    group_bonds_with_h_atomtype = get_atoms_bond(groups_bonds_with_h,atomtypedic)
                    groups_bonds_without_h_atomtype = get_atoms_bond(groups_bonds_without_h,atomtypedic)
                    groups_angles_with_h_atomtype = get_atoms_bond(groups_angles_with_h,atomtypedic)
                    groups_angles_without_h_atomtype = get_atoms_bond(groups_angles_without_h,atomtypedic)
                    
                    oldtop = parmed.load_file('FeCP2_0_dry_gromacs.top')

                    for key in bondinfos:
                        if set(key) not in list(group_bonds_with_h_atomtype.values()):
                            if set(key) not in list(groups_bonds_without_h_atomtype.values()):
                                pair = []
                                value = bondinfos[key]
                                for i in key:
                                    for j in atomtypedic:
                                        if i == atomtypedic[j]:
                                            pair.append(j)
                                print(pair)
                    
                                atomI = oldtop.atoms[pair[0]-1]
                                atomJ = oldtop.atoms[pair[1]-1]
                                fc = value[0]
                                eq = value[1]
                                print(fc,eq)
                                
                                new_type = parmed.BondType(fc,eq)
                                new_bond= parmed.Bond(atomI,atomJ,type=new_type)
                                oldtop.bonds.insert(0,new_bond)
                    print(oldtop.bonds)
                    oldtop.save(newprmtop)



a = checkff('FeCP2_0')
pdbinfos,lignameinfos= a.readpdb()
newpdbinfos,atomtypedic = a.readmol2(lignameinfos,pdbinfos)
#print(newpdbinfos,atomtypedic)
bondinfos, angleinfos, dihedralinfos = a.readfrcmod()
#print(bondinfos,angleinfos)
a.gen_new_prmtop(oldprmtop='FeCP2_0_dry.prmtop',newprmtop='FeCP2_0_dry_new4.prmtop',missingbond='missingbonds.txt',
                 atomtypedic=atomtypedic,bondinfos=bondinfos,angleinfos=angleinfos)



                    

                    
                   



                
        




        







