from openbabel import openbabel
import subprocess
import numpy as np
from ase.io import read, write
from molSimplify.Classes.mol3D import mol3D
from molSimplify.Classes.ligand import ligand_breakdown

import os
import re
import sys
import getopt
from glob import glob
from rdkit import Chem

atom_zoo = ['H','C','N','O','CL','BR','I','P','S']

class AtomInfo():
    def __init__(self,atomnum,atomtype,coordinate):
        self.atomnum = atomnum
        self.atomtype = atomtype
        self.coordinate = coordinate
"""
        class AtomInfo:
            def __init__(self, atomnumber, coordinate, atomtype):
                self.atomnumber = atomnumber
                self.coordinate = coordinate
                self.atomtype = atomtype
"""
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
    def __init__(self, xyzfile='', filename = None, metal_charge=None, metal_ID='',chargefile=None,mode='A',spinmult=None): #liglist='', denticity='',ligcons='',atomsinfo='',ligand_charge='',metal_name=''):
        self.xyzfile = xyzfile
        self.metal_charge = metal_charge
        self.metal_ID = metal_ID
        self.filename = filename
        self.xyzfile = filename + '.xyz'
        self.chargefile = chargefile
        self.mode = mode
        self.spinmult = spinmult
    
    def inputCheck(self):
        if self.filename == None:
            print('Error: No xyz file of ligand-metal complex is assigned')
            sys.exit()
        elif os.path.exists(self.filename  + '.xyz') == False:
            print('Error: can not find ' + self.filename  + '.xyz ' )
            sys.exit()
        if self.metal_charge == None:
            print('Error: No charge is assigned to metal')
            sys.exit()
        if self.mode not in ['A','M']:
            print('Error: invalid mode for ligand charge assignment')

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
        metal_list = complex_mol.findMetal()
        if len(metal_list) > 1:
            print("Error: AutoMCPB can't deal with multi-metal systems")
            sys.exit()
        elif len(metal_list) == 0:
            print('Error: AutoMCPB There is no metal in the xyz file, please check the xyz input file.')
        elif len(metal_list) == 1:
            self.metal_ID = metal_list[0]
            metal_type = self.atomsinfo[self.metal_ID].atomtype
            print('Find Metal ID:' + str(self.metal_ID) + ' ' + metal_type)
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
                print('Error: AutoMCPB does not support the atom ' + atom.atomtype)
                check = 'No'
                sys.exit()
        if check == 'Yes':
            print('All atom types are supported by AutoMCPB')
        
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
        ligands_breakdown_infos = ligand_breakdown(complex_mol)
        liglist = ligands_breakdown_infos[0]
       # denticity = ligands_breakdown_infos[1]
        ligcons =  ligands_breakdown_infos[2]
        self.liglist = liglist
        self.ligcons = ligcons
        #self.denticity = denticity
        if len(liglist) == 0:
            print('Error: No ligand is identified')
            sys.exit()
        else:
            print('In', self.xyzfile, ':')
            for ID, ligand_atoms in enumerate(liglist):
                print('Ligand ' + str(ID), ligand_atoms, 'is identified' )
                if len(ligand_atoms) < 2:
                    print('Error:','Ligand ' + str(ID) + ' has too few atoms to generate mol2 in GAFF.')
                  #  sys.exit()
            print('atom :',ligcons, 'linked with the metal',self.metal_ID )
    
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
        extra_charge = 'N'
        
            
        sdf = ligandxyz.split('.xyz')[0]
        ofile = open(sdf + '_' + self.metal_name + '.xyz','w')
        cmd = 'obabel ' + ligandxyz + ' -O ' + sdf + '.sdf'
        with open(sdf + '_sdf.log', 'w') as f:
            subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
    
        heavyatom_of_xyz = []

        ### get heavy atom of xyz ###
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
            metal_name = self.metal_name[0] + str.lower(self.metal_name[1])
            metal_x = self.atomsinfo[self.metal_ID].coordinate[0]
            metal_y = self.atomsinfo[self.metal_ID].coordinate[1]
            metal_z = self.atomsinfo[self.metal_ID].coordinate[2]
            ofile.write(metal_name+ ' ' + str(metal_x) + ' ' + str(metal_y) + ' ' + str(metal_z))
        ofile.close()

        suppl = Chem.SDMolSupplier(sdf+'.sdf')
        mol = [mol for mol in suppl][0]
        print('In ',ligandxyz.split('.xyz')[0],':')

        total_radical_number = 0
        for newID,atom in enumerate(mol.GetAtoms()):
            radical = atom.GetNumRadicalElectrons()
         #   print('atom', heavyatom_of_xyz[newID].atomnum,heavyatom_of_xyz[newID].atomtype, 'radical num is ', radical)
            if radical > 0:
                atoms_with_radical[heavyatom_of_xyz[newID].atomnum] = radical
                total_radical_number = total_radical_number + radical
        ligand_metal = mol3D()
        ligand_metal.readfromxyz(sdf + '_' + self.metal_name + '.xyz')
        ligand_metal_breakdown_infos = ligand_breakdown(ligand_metal)

        for atom in ligand_metal_breakdown_infos[2]:
            for ID in atoms_with_radical:
                if ID in atom:
                    charge_of_ligand = charge_of_ligand + atoms_with_radical[ID]*(-1)
                    print('atom',ID, 'withdraw', str(atoms_with_radical[ID]),'electron from metal')
                else:
                    extra_charge = 'Y'
                    print('atom',ID, 'have extra charges, please set the charge of each ligand in charge.in')
       
        if total_radical_number == 0:
            charge_of_ligand = 0
        print(sdf,'charge is',charge_of_ligand )
        return charge_of_ligand,extra_charge

    def ligand_breakdown(self):
        r"""
        breakdown the ligands from 
        parameters
        ----------
        None

        Returns             
        -------
        self.ligand_charge_dic, dictionary, key: ligand name like: 'LG1', 'LG2', self.ligand_charge_dic['LG1'] = charge, int
        self.extracharges, list of existence extracharge of each ligands, eg. ['N','Y','N' ..]
        """
      #  ligand_charge_inputfile = 'charge.in'
        ligand_charge_dic = {}
        extra_charges = []
      #  ofile_charge = open(ligand_charge_inputfile,'w')
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
            
            charge,extra_charge = self.check_radical(ligandname +'.xyz')
            extra_charges.append(extra_charge)
            ligand_charge_dic[ligandname] = charge
          #  ofile_charge.write(ligandname + ' ' +  str(charge) + '\n')          

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

        self.extracharges = extra_charges
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
        cmd = '$AMBERHOME/bin/'+'metalpdb2mol2.py -i ' + metalname + '.pdb'+ ' -o ' + metalname + '.mol2' + ' -c ' + str(self.metal_charge)
        subprocess.call(cmd, shell=True)
    
    def generate_mol2_by_auto(self):
        r"""
        automatically generate the mol2 file if self.mode = 'A'
        
        Parameters
        ----------
        None   
 
        Returns
        -------
        None
        """
        for ID, ligand in enumerate(self.liglist):
            ligandname = 'LG' + str(ID)
            ligand_charge_dic = self.ligand_charge_dic
            charge = ligand_charge_dic[ligandname]
            if len(ligand) > 1:
                cmd = '$AMBERHOME/bin/'+'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c bcc -pf y -nc ' + str(charge) + ' -m 2'
                print('***Generating the', ligandname + '.mol2','file...')
                with open(ligandname + '_antechamber_generate_mol2.log', 'w') as f:
                    subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
            elif len(ligand) == 1:
                cmd = '$AMBERHOME/bin/'+'metalpdb2mol2.py -i ' + ligandname + '___.pdb'+ ' -o ' + ligandname + '.mol2' + ' -c ' + str(charge)
                print('***Generating the', ligandname + '.mol2','file...')
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
                cmd = '$AMBERHOME/bin/'+'antechamber -fi pdb -fo mol2 -i ' + ligandname +'___.pdb' + ' -o ' + ligandname + '.mol2'  + ' -c bcc -pf y -nc ' + str(charge) + ' -m 2'
                print('***Generating the', ligandname + '.mol2','file...')
                with open(ligandname + '_antechamber_generate_mol2.log', 'w') as f:
                    subprocess.call(cmd, shell=True, stdout=f, stderr=subprocess.STDOUT)
            
            elif length == 1:
                cmd = '$AMBERHOME/bin/'+'metalpdb2mol2.py -i ' + ligandname + '___.pdb'+ ' -o ' + ligandname + '.mol2' + ' -c ' + str(charge)
                print('***Generating the', ligandname + '.mol2','file...')
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
                cmd = '$AMBERHOME/bin/'+ 'parmchk2 -i ' + ligandname + '.mol2' + ' -o ' + ligandname +'.frcmod' +' -f mol2'
                subprocess.call(cmd,shell=True)
                self.modifiy_pdb(ligandname +'___.pdb',ligandname + '.mol2',ligandname  +'_temp.pdb')
                cmd = 'cat '+ ligandname  + '_temp.pdb  | grep ' + ligandname + ' > ' + ligandname  + '.pdb'  
                subprocess.call(cmd,shell=True)   
                
            else:
                print(ligandname + '.mol2', 'file', 'is note generated by antechamber, please check', ligandname + '_antechamber_generate_mol2.log')

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
        subprocess.call(cmd,shell=True)
        cmd = '$AMBERHOME/bin/'+ 'pdb4amber -i ' + temppdb + ' -o ' +  self.filename + '_final.pdb'
        with open(ligandname + '_pdb4amber.log', 'w') as f:
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
        pdbin = self.filename + '_final.pdb'
        xyzout = self.filename + '_final.xyz'
        write(xyzout, read(pdbin))        
        new_complex_mol = mol3D()
        new_complex_mol.readfromxyz(self.filename + '_final.xyz')
        metalID = new_complex_mol.findMetal()
        self.new_complex_mol = new_complex_mol
        pairs = ligand_breakdown(new_complex_mol)[2]

        add_bonded_pairs = 'add_bonded_pairs '
        for denticitys in pairs:
            for denticity in denticitys:
                metal_denticity = str(metalID[0]+1) + '-' + str(denticity+1) + ' '
                add_bonded_pairs = add_bonded_pairs + metal_denticity
        
        self.add_bonded_pairs = add_bonded_pairs
        
    def generate_MCPB_input(self):
        metalname = self.atomsinfo[self.metal_ID].atomtype
        metalID = self.new_complex_mol.findMetal()
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
            f.write(original_pdb + '\n')
            f.write('group_name ' + self.filename + '\n')
            f.write('cutoff 2.8\n')
            f.write(ion_ids + '\n')
            f.write(self.add_bonded_pairs + '\n')
            f.write(ion_mol2files+ '\n')
            f.write(naa_mol2files+ '\n')
            f.write(frcmod_files + '\n')
            f.write('force_field  ff14SB\n')
            f.write('gaff 1\n')
            if self.spinmult == None:
                print('the spin is default by the number of electrons')
            else:
                if self.spinmult is isinstance(self.spinmult, int):
                    print('the spin is set as '+ str(self.spinmult))
                    f.write('lgmodel_spin ' + str(self.spinmult)+'\n')
                    f.write('smmodel_spin ' + str(self.spinmult)+'\n')
                else:
                    print('Error: The input for spins is not int')
    
    def get_connection_info(self,)
    
    def build(self):
        self.inputCheck()
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
     #   print("The autoMCPB script has been finished successfully")
    
def startautoMCPB(argumentList):
    options = "hn:c:u:m:f:"
    long_options = ["help",'filename=','metal_charge=','spin=','mode=','chargefile=']
    arguments, values = getopt.getopt(argumentList,options,long_options)
    filename = None
    metal_charge = None
    mode='A'
    spinmult=None
    chargefile = None
    print(arguments)
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
                                        FG1 -1 # ligand name charge 
            Example:
                python autoMCPB-smiles-v3.py -n 283-2-water --mode A
                '''
            print(message)
            sys.exit()
        elif currentArgument in ("-n","--filename"):
            print ("metal_complex filename", currentValue)
            filename = str(currentValue)
        elif currentArgument in ("-c",'--metal_charge'):
            metal_charge = int(currentValue)
        elif currentArgument in ('-u','--spin'):
            spinmult=int(currentValue)
        elif currentArgument in ('-m','--mode'):
            mode = str(currentValue)
            if mode not in ['A','M']:
                print('Error: invalid argument for --mode')
        elif currentArgument in ('-f','--chargefile'):
            chargefile = str(currentValue)
            print('read the charge from', chargefile)

    builder = AutoMCPB(filename=filename,metal_charge=metal_charge, spinmult=spinmult,
                       mode=mode,chargefile=chargefile)

    builder.build()

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startautoMCPB(argumentList)

#print(a.ligcons)
