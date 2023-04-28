#--------------------------------------------------------------------------------------------------#
# A tool class for frcmod file generation, modification and analysis
# author: Fangning Ren (2023-04-27)
# path: autosolvate/utils/frcmod.py
#--------------------------------------------------------------------------------------------------#

import os
import re
import shutil
import subprocess
from typing import List, Dict, Tuple, TextIO, Optional

import numpy as np

# frcmod data are stored in the following format:
# format of self.data:
# {
#     "MASS": {
#         "C": [12.011,],
#         "H": [1.008,],
#         ...
#     },
#     "BOND": {
#         ("C", "C"): [300.0, 1.54],
#         ("C", "H"): [300.0, 1.09],
#         ...
#     },
#     "ANGL": {
#         ("C", "C", "C"): [50.0, 120.0],
#         ("C", "C", "H"): [50.0, 109.5],
#         ...
#     },
#     "DIHE": {
#         ("C", "C", "C", "C"): [1, 1.2, 180.0, 2.0],
#         ("C", "C", "C", "H"): [1, 1.2, 180.0, 2.0],
#         ...
#     },
#     "IMPROPER": {
#         ("C", "C", "C", "C"): [1.1, 180.0, 2.0],
#         ("C", "C", "C", "H"): [1.1, 180.0, 2.0],
#         ...
#     },
#     "NONB": {
#         "C": [1.9080, 0.0860],
#         "H": [1.4870, 0.0157],
#         ...
#     },
#     "HBON": {
#         ("C", "C"): [0.0, 0.0],
#         ("C", "H"): [0.0, 0.0],
#         ...
#     },


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


if __name__ == "__main__":
    frcmod = FrcmodFile("ppi.frcmod")
    urcmod = FrcmodFile("ppi copy.frcmod")
    frcmod.update(urcmod)
    frcmod.write("ppi-updated.frcmod")