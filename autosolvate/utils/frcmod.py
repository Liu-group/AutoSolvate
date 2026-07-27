#--------------------------------------------------------------------------------------------------#
# A tool class for frcmod file generation, modification and analysis
# author: Fangning Ren (2024-04-16)
# path: autosolvate/utils/frcmod.py
#--------------------------------------------------------------------------------------------------#

import os
import re
import shutil
import subprocess
from typing import List, Dict, Tuple, TextIO, Optional
from collections import OrderedDict
from copy import deepcopy


import numpy as np

from parmed.parameters import ParameterSet
from parmed.amber import AmberParameterSet, AmberParm
from parmed.gromacs import GromacsTopologyFile
from parmed.exceptions import FormatNotFound

"""
frcmod data are stored in the following format:
format of self.data:
{
    "MASS": {
        "C": [12.011,],
        "H": [1.008,],
        ...
    },
    "BOND": {
        ("C", "C"): [300.0, 1.54],
        ("C", "H"): [300.0, 1.09],
        ...
    },
    "ANGL": {
        ("C", "C", "C"): [50.0, 120.0],
        ("C", "C", "H"): [50.0, 109.5],
        ...
    },
    "DIHE": {
        ("C", "C", "C", "C"): [1, 1.2, 180.0, 2.0],
        ("C", "C", "C", "H"): [1, 1.2, 180.0, 2.0],
        ...
    },
    "IMPROPER": {
        ("C", "C", "C", "C"): [1.1, 180.0, 2.0],
        ("C", "C", "C", "H"): [1.1, 180.0, 2.0],
        ...
    },
    "NONB": {
        "C": [1.9080, 0.0860],
        "H": [1.4870, 0.0157],
        ...
    },
    "HBON": {
        ("C", "C"): [0.0, 0.0],
        ("C", "H"): [0.0, 0.0],
        ...
    },
"""

class FrcmodFile:

    def __init__(self, frcmod:str):
        self.frcmod = frcmod
        self.lines = []
        self.data = {
            "MASS"      : {},
            "BOND"      : {},
            "ANGL"      : {},
            "DIHE"      : {},
            "IMPROPER"  : {},
            "NONBON"    : {},
            "HBOND"     : {},
        }
        self.read()
        self.parse()

    def read(self):
        with open(self.frcmod, "r") as f:
            self.lines = f.readlines()
   
    def parse(self):
        curlineindex = 0
        for i, line in enumerate(self.lines):
            if line.startswith("MASS"):
                curlineindex += i + 1
                break

        for i, line in enumerate(self.lines[curlineindex:]):
            if line.startswith("BOND"):
                curlineindex += i + 1
                break
            line = line.strip()
            if not line or len(line) <= 2:
                continue
            line = line[:2].replace("-", " ") + line[2:]
            words = line.split()
            if len(words) < 2:
                continue
            self.data["MASS"][words[0]] = [float(words[1]),]

        for i, line in enumerate(self.lines[curlineindex:]):
            if line.startswith("ANGL") or line.startswith("ANGLE"):
                curlineindex += i + 1
                break
            line = line.strip()
            if not line or len(line) <= 5:
                continue
            line = line[:5].replace("-", " ") + line[5:]
            words = line.split()
            if len(words) < 4:
                continue
            bkey = (words[0], words[1]) 
            self.data["BOND"][bkey] = [float(words[2]), float(words[3])]

        for i, line in enumerate(self.lines[curlineindex:]):
            if line.startswith("DIHE") or line.startswith("DIHEDRAL"):
                curlineindex += i + 1
                break
            line = line.strip()
            if not line or len(line) <= 8:
                continue
            line = line[:8].replace("-", " ") + line[8:]
            words = line.split()
            if len(words) < 5:
                continue
            akey = (words[0], words[1], words[2])
            self.data["ANGL"][akey] = [float(words[3]), float(words[4])]

        for i, line in enumerate(self.lines[curlineindex:]):
            if line.startswith("IMPR") or line.startswith("IMPROPER"):
                curlineindex += i + 1
                break
            line = line.strip()
            if not line or len(line) <= 11:
                continue
            line = line[:11].replace("-", " ") + line[11:]
            words = line.split()
            if len(words) < 7:
                continue
            dkey = (words[0], words[1], words[2], words[3])
            self.data["DIHE"][dkey] = [int(words[4]), float(words[5]), float(words[6]), float(words[7])]

        for i, line in enumerate(self.lines[curlineindex:]):
            if line.startswith("NONB") or line.startswith("NONBON"):
                curlineindex += i + 1
                break
            line = line.strip()
            if not line or len(line) <= 11:
                continue
            line = line[:11].replace("-", " ") + line[11:]
            words = line.split()
            if len(words) < 4:
                continue
            ikey = (words[0], words[1], words[2], words[3])
            self.data["IMPROPER"][ikey] = [float(words[4]), float(words[5]), float(words[6])]

        for i, line in enumerate(self.lines[curlineindex:]):
            if line.startswith("HBON") or line.startswith("HBOND"):
                curlineindex += i + 1
                break
            line = line.strip()
            if not line or len(line) <= 2:
                continue
            line = line[:2].replace("-", " ") + line[2:]
            words = line.split()
            if len(words) < 3:
                continue
            nkey = words[0]
            self.data["NONBON"][nkey] = [float(words[1]), float(words[2])]

        for i, line in enumerate(self.lines[curlineindex:]):
            line = line.strip()
            if not line or len(line) <= 2:
                continue
            line = line[:2].replace("-", " ") + line[2:]
            words = line.split()
            if len(words) < 4:
                continue
            hkey = (words[0], words[1])
            self.data["HBOND"][hkey] = [float(words[2]), float(words[3])]

    def update(self, frcmod):
        for key in frcmod.data:
            self.data[key].update(frcmod.data[key])

    def get_mdgx_labels(self):
        bonds = [("{:<2s}".format(key[0]), "{:<2s}".format(key[1])) for key in self.data["BOND"].keys()]
        angles = [("{:<2s}".format(key[0]), "{:<2s}".format(key[1]), "{:<2s}".format(key[2])) for key in self.data["ANGL"].keys()]
        dihedrals = [("{:<2s}".format(key[0]), "{:<2s}".format(key[1]), "{:<2s}".format(key[2]), "{:<2s}".format(key[3])) for key in self.data["DIHE"].keys()]
        impropers = [("{:<2s}".format(key[0]), "{:<2s}".format(key[1]), "{:<2s}".format(key[2]), "{:<2s}".format(key[3])) for key in self.data["IMPROPER"].keys()]
        return bonds, angles, dihedrals, impropers

    def write(self, frcmod:str):
        f = open(frcmod, "w")
        for key, value in self.data.items():
            f.write(f"\n{key}\n")
            for k, v in value.items():
                term = "-".join(["{:<2s}"] * len(k)).format(*k)
                if key == "MASS":
                    f.write(f"{term} {v[0]:>10.4f}\n")
                elif key == "BOND":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f}\n")
                elif key == "ANGL" or key == "ANGLE":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f}\n")
                elif key == "DIHE" or key == "DIHEDRAL":
                    f.write(f"{term} {v[0]:>10d} {v[1]:>10.4f} {v[2]:>10.4f} {v[3]:>10.4f}\n")
                elif key == "IMPROPER" or key == "IMPROPER":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f} {v[2]:>10.4f}\n")
                elif key == "NONB" or key == "NONBON":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f}\n")
                elif key == "HBON" or key == "HBOND":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f}\n")
            f.write("\n")
        f.close()

    def write_with_prm_label(self, frcmod:str, 
            fit_b_length:bool=True , 
            fit_b_fconst:bool=True ,
            fit_a_angle: bool=True ,
            fit_a_fconst:bool=True ,
            fit_d_phase: bool=True ,
            fit_d_fconst:bool=True ,
            fit_i_fconst:bool=False,
            fit_i_phase: bool=False,
            fit_n_radius:bool=False,
            fit_n_well:  bool=False,
            fit_h_radius:bool=False,
            fit_h_well:  bool=False):
        bondlabel = "PRM"
        if fit_b_fconst:
            bondlabel += " 1"
        if fit_b_length:
            bondlabel += " 2"
        if not fit_b_fconst and not fit_b_length:
            bondlabel = ""
        anglelabel = "PRM"
        if fit_a_fconst:
            anglelabel += " 1"
        if fit_a_angle:
            anglelabel += " 2"
        if not fit_a_fconst and not fit_a_angle:
            anglelabel = ""
        dihlabel = "PRM"
        if fit_d_fconst:
            dihlabel += " 2"
        if fit_d_phase:
            dihlabel += " 3"
        if not fit_d_fconst and not fit_d_phase:
            dihlabel = ""
        imlabel = "PRM"
        if fit_i_fconst:
            imlabel += " 1"
        if fit_i_phase:
            imlabel += " 2"
        if not fit_i_fconst and not fit_i_phase:
            imlabel = ""
        nonlabel = "PRM"
        if fit_n_radius:
            nonlabel += " 1"
        if fit_n_well:
            nonlabel += " 2"
        if not fit_n_radius and not fit_n_well:
            nonlabel = ""
        hbonlabel = "PRM"
        if fit_h_radius:
            hbonlabel += " 1"
        if fit_h_well:
            hbonlabel += " 2"
        if not fit_h_radius and not fit_h_well:
            hbonlabel = ""

        f = open(frcmod, "w")
        for key, value in self.data.items():
            f.write(f"\n{key}\n")
            for k, v in value.items():
                term = "-".join(["{:<2s}"] * len(k)).format(*k)
                if key == "MASS":
                    f.write(f"{term} {v[0]:>10.4f}\n")
                elif key == "BOND":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f} {bondlabel}\n")
                elif key == "ANGL" or key == "ANGLE":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f} {anglelabel}\n")
                elif key == "DIHE" or key == "DIHEDRAL":
                    f.write(f"{term} {v[0]:>10d} {v[1]:>10.4f} {v[2]:>10.4f} {v[3]:>10.4f} {dihlabel}\n")
                elif key == "IMPROPER" or key == "IMPROPER":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f} {v[2]:>10.4f} {imlabel}\n")
                elif key == "NONB" or key == "NONBON":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f} {nonlabel}\n")
                elif key == "HBON" or key == "HBOND":
                    f.write(f"{term} {v[0]:>10.4f} {v[1]:>10.4f} {hbonlabel}\n")
            f.write("\n")
        f.close()

# functions utilizing parmed to modify the frcmod file
def remove_zero_dihedrals(frcmod:str, out:str=None):
    """
    Remove dihedral terms with per(periodicity) = 0.0, so programs as OpenMM will not complain about it
        for dihedral with per = 0.0 and phi_k != 0.0, remove this term
        for dihedral with per = 0.0 and phi_k = 0.0, change the per to 1.0 (this will not affect the energy calculation)
    """
    if not out:
        out = frcmod
        shutil.copy(frcmod, os.path.splitext(out)[0] + "-original.frcmod")
    a = AmberParameterSet(frcmod)
    for key, dihlist in a.dihedral_types.items():
        id_of_dih_to_be_removed = []
        for dih in dihlist:
            if dih.phi_k == 0.0 and dih.per == 0.0:
                dih.per = 1.0
            elif dih.phi_k > 0.0 and dih.per == 0.0:
                id_of_dih_to_be_removed.append(id(dih))
        for dih_id in id_of_dih_to_be_removed:
            j = 0
            while j < len(dihlist):
                if id(dihlist[j]) == dih_id:
                    del dihlist[j]
                else:
                    j += 1
    a.write(out)


# functions for pre-processing gromacs itp file
def read_gromacs_itp(itp_path:str):
    """
    Read gromacs itp file and return the atom types and charges
    """
    current_section = "" 
    section_contents = OrderedDict()
    with open(itp_path, 'r') as f:
        for line in f: 
            line = line.strip()
            if line.startswith(';') or line == '':
                continue
            line = line.split(";")[0]
            if line.startswith('['):
                line = line.replace('[', '').replace(']', '')
                current_section = line.strip()
                section_contents[current_section] = []
                continue
            else:
                section_contents[current_section].append(line)
    return section_contents

def parse_gromacs_itp(section_contents):
    molecule_info = []      # {"resname":str, "nrexcl":int}
    atom_info = []
    bond_info = []
    angle_info = []
    dihedral_info = []
    for line in section_contents["moleculetype"]:
        data = line.split()
        if len(data) <= 1:
            continue
        resname = data[0]
        nrexcl = int(data[1])
        molecule_info.append({"resname":resname, "nrexcl":nrexcl})
    for line in section_contents["atoms"]:
        data = line.split()
        if len(data) <= 8:
            continue
        atom_info.append({
            "nr":int(data[0]), "type":data[1], "resnr":int(data[2]), "residue":data[3], 
            "atom":data[4], "cgnr":int(data[5]), "charge":float(data[6]), "mass":float(data[7])}
        )
    for line in section_contents["bonds"]:
        data = line.split()
        if len(data) <= 5:
            continue
        bond_info.append({
            "ai":int(data[0]), "aj":int(data[1]), "funct":int(data[2]), 
            "params":list(map(float, data[3:]))
        })
        if int(data[2]) != 1:
            print("Warning: bond function type between atom {} and atom {} is not 1".format(data[0], data[1]))
    for line in section_contents["angles"]:
        data = line.split()
        if len(data) <= 6:
            continue
        angle_info.append({
            "ai":int(data[0]), "aj":int(data[1]), "ak":int(data[2]), "funct":int(data[3]), 
            "params":list(map(float, data[4:]))
        })
        if int(data[3]) != 1:
            print("Warning: angle function type between atom {} and atom {} and atom {} is not 1".format(data[0], data[1], data[2]))
    for line in section_contents["dihedrals"]:
        data = line.split()
        if len(data) <= 7:
            continue
        dihedral_info.append({
            "ai":int(data[0]), "aj":int(data[1]), "ak":int(data[2]), "al":int(data[3]), "funct":int(data[4]), 
            "params":list(map(float, data[5:]))
        })
        if int(data[4]) not in [1,2,3,4,5,9]:
            print("Warning: dihedral function type between atom {} and atom {} and atom {} and atom {} is not 1, 2, 3, 4, 5 or 9".format(data[0], data[1], data[2], data[3]))
    
    return molecule_info, atom_info, bond_info, angle_info, dihedral_info

def convert_dihedral_5to3(dihedral_list):
    type_5_to_3 = np.array([
        [ 0.5,  1.0,  0.5,  0.0], 
        [-0.5,  0.0,  1.5,  0.0],
        [ 0.0, -1.0,  0.0,  4.0],
        [ 0.0,  0.0, -2.0,  0.0],
        [ 0.0,  0.0,  0.0, -4.0],
        [ 0.0,  0.0,  0.0,  0.0],
    ])
    for dihedral in dihedral_list:
        funct_type = dihedral['funct']
        if funct_type != 5:
            continue
        params = dihedral['params']
        params = np.array(params)
        if len(params) < 4:
            params = np.concatenate([params, np.zeros(4-len(params))])
        elif len(params) > 4:
            params = params[:4]
            print("Warning: Fourier series coefficients for dihedral {} {} {} {} has more than 4 parameters. Truncated into 4.".format(dihedral['ai'], dihedral['aj'], dihedral['ak'], dihedral['al']))
        params = np.dot(type_5_to_3, params)
        params = np.round(params, 4)
        dihedral['funct'] = 3
        dihedral['params'] = params.tolist()
    return dihedral_list

def reformat_lines(section_content, target_section, new_data):
    if target_section not in ["atoms", "bonds", "angles", "dihedrals"]:
        return section_content
    new_lines = []
    if target_section == "atoms":
        for atom in new_data:
            new_lines.append("{nr:>5} {type:<10} {resnr:>5} {residue:<10} {atom:<10} {cgnr:>5} {charge:>8.4f} {mass:>8.4f}".format(**atom))
    elif target_section == "bonds":
        for bond in new_data:
            new_lines.append("{ai:>5} {aj:>5} {funct:>5}".format(**bond) + " ".join(["{:>10.4f}".format(p) for p in bond['params']]))
    elif target_section == "angles":
        for angle in new_data:
            new_lines.append("{ai:>5} {aj:>5} {ak:>5} {funct:>5} ".format(**angle) + " ".join(["{:>10.4f}".format(p) for p in angle['params']]))
    elif target_section == "dihedrals":
        for dihedral in new_data:
            new_lines.append("{ai:>5} {aj:>5} {ak:>5} {al:>5} {funct:>5} ".format(**dihedral) + " ".join(["{:>10.4f}".format(p) for p in dihedral['params']]))
    section_content[target_section] = new_lines
    return section_content

def write_gromacs_itp(itp_path:str, section_contents:OrderedDict):
    with open(itp_path, 'w') as f:
        for section_name, section_data in section_contents.items():
            f.write("[  {}  ]\n".format(section_name))
            for line in section_data:
                f.write(line + "\n")
            f.write("\n")

def convert_dihedral_5to3_main(itp_path:str):
    """
    convert the dihedral type 5 into 3 in the gromacs itp file
    """
    shutil.copyfile(itp_path, os.path.splitext(itp_path)[0] + "-type5dih.itp")
    section_contents = read_gromacs_itp(itp_path)
    molecule_info, atom_info, bond_info, angle_info, dihedral_info = parse_gromacs_itp(section_contents)
    dihedral_info = convert_dihedral_5to3(dihedral_info)
    section_contents = reformat_lines(section_contents, "dihedrals", dihedral_info)
    write_gromacs_itp(itp_path, section_contents)

def read_molecule_name(itp_file_path):
    with open(itp_file_path, 'r') as file:
        lines = file.readlines()
    lines = [line.strip() for line in lines]
    lines = [line for line in lines if not line.startswith(';')]
    moleculesectionindex = -1
    for i, line in enumerate(lines):
        if line.startswith("["):
            section = line.replace("[", "").replace("]", "").strip()
            if section == "moleculetype":
                moleculesectionindex = i
                break
    if moleculesectionindex == -1:
        return "UNK"
    molecule_name = lines[moleculesectionindex + 1].split()[0]
    return molecule_name

def write_single_molecule_top(itp_pathes:List[str], output_path:str):
    """
    create a gromacs top file that only contains single molecule so parmed can read it 
    """
    mol_res = "UNK"
    for itp_path in itp_pathes:
        if not os.path.exists(itp_path):
            raise FileNotFoundError("File not found: %s" % itp_path)
        molecule_res = read_molecule_name(itp_path)
        if molecule_res != "UNK":
            mol_res = molecule_res
            break
    with open(output_path, 'w') as file:
        file.write("#define _FF_OPLS\n")
        file.write("#define _FF_OPLSAA\n")
        file.write("[ defaults ]\n")
        file.write(";nbfunc  comb-rule   gen-pairs   fudgeLJ  fudgeQQ\n")
        file.write("1         3          yes        0.5      0.5\n")
        file.write("\n")
        for itp_path in itp_pathes:
            file.write("#include \"%s\"\n" % itp_path)
        file.write("\n")
        file.write("[ system ]\n")
        file.write("Neat %s\n" % mol_res)
        file.write("\n")
        file.write("[ molecules ]\n")
        file.write("%s 1\n" % mol_res)

def itp_to_frcmod_mol2(itp_path:str, frcmod_path:str, mol2_path:str, xyz_path:str=None):
    """
    convert gromacs itp file to frcmod file and a mol2 file
    """
    itp_pathes = [itp_path,] if isinstance(itp_path, str) else itp_path
    for itp_path in itp_pathes:
        convert_dihedral_5to3_main(itp_path)
    toppath = os.path.splitext(itp_path[0])[0] + ".top"
    write_single_molecule_top(itp_pathes, toppath)
    try:
        gromacs_top = GromacsTopologyFile(toppath, parametrize=True, xyz=xyz_path)
    except FormatNotFound:
        gromacs_top = GromacsTopologyFile(toppath, parametrize=True)
    gromacs_top.save("tmp.prmtop", format="amber", overwrite=True)
    gromacs_top.save(mol2_path, format="mol2", overwrite=True)
    amber_parms = AmberParameterSet.from_structure(gromacs_top)
    amber_parms.write(frcmod_path)
    remove_zero_dihedrals(frcmod_path)


if __name__ == "__main__":
    frcmod = FrcmodFile("ppi.frcmod")
    urcmod = FrcmodFile("ppi copy.frcmod")
    frcmod.update(urcmod)
    frcmod.write("ppi-updated.frcmod")