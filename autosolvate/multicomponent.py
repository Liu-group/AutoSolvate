#--------------------------------------------------------------------------------------------------#
# multicomponent.py. 
# Description: 
#   This module can handle structure files containing multiple molecules.
#   1. Generates the lib and frcmod file for each separate molecule when it is not h2o or amino acids.
#   2. Generates the lib, prmtop and inpcrd files of this whole structure.
# Update 2022-02-04:
#   1. Can generate lib, mol2 and prmtop file for an xyz with mutiple fragments with charges
#   2. Can create solvent box for xyz file with multiple fragments.
#   3. Able to rearrange shuffled pdb files into ordered form that antechamber can process.
# author: Fangning Ren (2022-02-04) 
# path: autosolvate/multicomponent.py
#--------------------------------------------------------------------------------------------------#
import getopt, sys, os
import subprocess
from openbabel import pybel
from openbabel import openbabel as ob

from autosolvate.autosolvate import *

def check_multicomponent(filename:str):
    mol = pybel.readfile(os.path.splitext(filename)[-1][1:], filename).__next__()
    mol_obmol = mol.OBMol
    fragments = mol_obmol.Separate()
    return len(fragments) > 1

def splitpdb(pdbname:str):
    with open(pdbname, "r") as f:
        lines = f.readlines()
    atomstart = 0
    header = []
    for i, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            atomstart = i
            break
        header.append(line)
    header = "".join(header)
    atomlines = []
    atomend = atomstart
    for i in range(atomstart, len(lines)):
        line = lines[i]
        if line.startswith("TER"):
            continue
        if not line.startswith("ATOM") and not line.startswith("HETATM"):
            atomend = i
            break
        atomid = int(line[6:12])-1
        # print(i, atomid, len(atomlines))
        assert len(atomlines) == atomid
        atomlines.append(line)
    middles = []
    conectstart = atomend
    for i in range(atomend, len(lines)):
        if lines[i].startswith("CONECT"):
            conectstart = i
            break
        middles.append(lines[i])
    conectend = conectstart
    bonds = []
    for i in range(conectstart, len(lines)):
        if not lines[i].startswith("CONECT"):
            conectend = i
            break
        bondline = tuple(map(int, lines[i].split()[1:]))
        bonds.append(bondline)
    tailer = []
    for i in range(conectend, len(lines)):
        tailer.append(lines[i])
    tailer = "".join(tailer)
    return header, atomlines, bonds, tailer

class MulticomponentParamsBuilder():
    def __init__(self, xyzfile: str, 
    name="", resname="", charge=0, spinmult=1, charge_method="resp", outputFile="", pre_optimize_fragments = False,
    srun_use=False, gaussianexe=None, gaussiandir=None, amberhome=None, deletefiles=False): 
        """
        Create amber parameter files for a single xyz or pdb file with multiple separate fragments

        Parameters
        ----------
        xyzfile : str
            structure file name, can be any structural files that openbabel recognizes.
        name : array_like, Optional. default: the base name of the provided structure file.
            Not used
        resname : array_like, Optional. default: Residue name provided in pdb file or assigned as UAA, UAB, UAC, etc.
            Residue names for each fragments. A list of strings of three capital letters. Its length should equal to the number of fragments in xyzfile. If this parameter is not given, the residues will be assigned by "U" plus "AB","AC",..."AZ", "BA"...
        charge : dict | array_like, Optional. default: 0
            Charge for each fragment. A list of integer with its length equal to the number of fragments, or a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as neutral. 
        spinmult : dict | array_like, Optional. default: 0
            Multiplicity for each fragment. A list of integer with its length equal to the number of fragments, or a dictionary with the three-letter name of the residue as the key and the corresponding charge as the value. If not given, all fragment will be considered as singlet. 
        outputFile : str, Optional, default='water_solvated'
            Filename-prefix for outputfiles
        pre_optimize_fragments : bool, Optional, default: False
            do geometry optimization with MMFF94 forcefield in OpebBabel before running antechamber
        srun_use : bool, Optional, default='False
            Run all commands with a srun prefix
        gaussianexe : str, Optional, default: g16
            name of the Gaussian executeble
        gaussiandir : str, Optional, default: $GAUSSIANDIR
            path of Gaussian
        amberhome : str, Optional, default: $AMBERHOME
            path of amber
        deletefiles : bool, Optional, default: False
            Delete all temporary files except the .prmtop and .inpcrd file of the pdb file provided.
        """
        self.names = name
        self.resnames = resname
        self.charges = charge
        self.netcharge = 0
        self.spinmults = spinmult
        self.charge_method = charge_method

        self.mainfilename = xyzfile
        self.basename, ext = os.path.splitext(self.mainfilename)
        filename = os.path.basename(self.mainfilename)
        
        # copy the input pdb to the current work dir
        if not os.path.exists(os.path.join(os.getcwd(), filename)):
            f = open(self.mainfilename, "r")
            content = f.read()
            f.close()
            self.mainfilename = os.path.join(os.getcwd(), filename)
            self.basename, ext = os.path.splitext(self.mainfilename)
            f = open(self.mainfilename, "w")
            f.write(content)
            f.close()
        if True:
            f = open(self.mainfilename, "r")
            content = f.read()
            f.close()
            f = open(self.basename + "-original" + ext, "w")
            f.write(content)
            f.close()

        self.pathname = os.path.dirname(self.basename)
        if self.pathname == "":
            self.pathname = os.getcwd()
        self.name = os.path.basename(self.basename)

        self.outputFile = outputFile
        self.srun_use = srun_use
        self.gaussianexe = gaussianexe
        self.gaussiandir = gaussiandir
        self.amberhome = amberhome
        self.pre_optimize_fragments = pre_optimize_fragments
        self.deletefiles = deletefiles

        self.aminoacidresidues = set(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"])

        self.getPDB()
        self.mol_obmol = pybel.readfile("pdb", self.mainfilename).__next__().OBMol
        self.fragnames = []
        self.fragmols = []
        self.newfragmols = []
        self.newfragpdbs = []
        self.newresiduenames = []

    def getPDB(self):
        """convert the input structure file to pdb"""
        basename, ext = os.path.splitext(self.mainfilename)
        ext = ext[1:]
        if ext != "pdb":
            subprocess.run(f"obabel -i {ext} {self.mainfilename} -o pdb > {basename}.pdb", shell = True)
        self.mainfilename = basename + ".pdb"

    def formatPDB(self, pdbfile:str):
        """format the pdb to the standard amber pdb using pdb4amber"""
        mainname, ext = os.path.splitext(pdbfile)
        os.rename(pdbfile, mainname + "-original.pdb")
        subprocess.run(f"pdb4amber -i {mainname}-original.pdb -o {pdbfile} -l {mainname}-format.log", shell = True)
        os.remove(f"{mainname}-original.pdb")

    def reorderPDB(self, pdbname:str):
        """
        rearrange order of atoms in pdb file. Atoms in each fragments after reshaping will have consecutive ids, with most of its original information preserved.
        Parameters
        ----------
        pdbname : str
            pdb file that need to be reordered
        """
        mol = pybel.readfile("pdb", pdbname).__next__()
        mol:pybel.Molecule
        mol_obmol = mol.OBMol
        mol_obmol:ob.OBMol
        n_atom = mol_obmol.NumAtoms()
        fragments = mol_obmol.Separate()
        newlabels = [-1 for i in range(n_atom)]
        baseindex = 0
        iscorrectorder = [True for i in range(len(fragments))]
        for i, fragment in enumerate(fragments):
            fragment:ob.OBMol
            n_atom_frag = fragment.NumAtoms()
            for j in range(n_atom_frag):
                a = fragment.GetAtom(j+1)
                newlabels[a.GetIdx()-1 + baseindex] = a.GetId()
                iscorrectorder[i] = (a.GetId() == a.GetIdx()-1 + baseindex)
            baseindex += n_atom_frag
        for flag in iscorrectorder:
            if not flag:
                print(f"WARNING: The atom order in PDB file {pdbname} is shuffled. Reordering...")
                break
        else:
            print("Atoms in each fragments have consecutive ids.")
        header, atomlines, bonds, tailer = splitpdb(pdbname)
        newlabelt = [-1 for i in range(len(newlabels))]

        basename = os.path.splitext(pdbname)[0]
        os.rename(pdbname, basename + "-not-ordered.pdb")
        with open(pdbname, "w") as f:
            f.write(header)
            for i, atomid in enumerate(newlabels):
                line = atomlines[atomid]
                f.write(line[0:6] + f"{i+1:>5d}" + line[11:])
            for i, lb in enumerate(newlabels):
                newlabelt[lb] = i
            for i in range(len(bonds)):
                bond = bonds[newlabels[i]]
                f.write("CONECT")
                for I in bond:
                    f.write(f"{newlabelt[I-1]+1:>5d}")
                f.write("\n")
            f.write(tailer)

    def checkParams(self):
        """check parameters, especially for the length of charge and multiplicity"""
        n_frag = len(self.fragmols)
        print(f"All fragments: {' '.join(self.fragnames)}")
        print(f"New fragments: {' '.join(self.newresiduenames)}")
        if isinstance(self.charges, int):
            print("WARNING: All charges are set to 0")
            self.charges = {frname:0 for frname in self.fragnames}
        if isinstance(self.spinmults, int):
            print("WARNING: All multiplicities are set to 1")
            self.spinmults = {frname:1 for frname in self.fragnames}
        if isinstance(self.charges, list):
            if len(self.charges) != len(self.fragnames):
                raise ValueError(f"Only charge for {len(self.charges)} fragments are provided. This file has {len(self.fragnames)} fragments!")
            self.charges = {frname:c for (frname, c) in zip(self.fragnames, self.charges)}
        if isinstance(self.spinmults, list):
            if len(self.spinmults) != len(self.fragnames):
                raise ValueError(f"Only multiplicity for {len(self.spinmults)} fragments are provided. This file has {len(self.fragnames)} fragments!")
            self.spinmults = {frname:c for (frname, c) in zip(self.fragnames, self.spinmults)}
        if isinstance(self.charges, dict) and isinstance(self.spinmults, dict):
            for frname in self.fragnames:
                if len(frname) > 3:
                    print(f"{frname} is a peptide with standard residues.")
                elif frname in self.newresiduenames:
                    if frname[0] == "U":
                        print(f"{frname} corresponds to a new fragment with auto-generated name.")
                    else:
                        print(f"{frname} corresponds to a new fragment.")
                    if frname not in self.charges:
                        print(f"WARNING: charge for {frname} not defined! Set it to 0 by default.")
                        self.charges[frname] = 0
                    else:
                        print(f"Set charge for {frname} to {self.charges[frname]}")
                    if frname not in self.spinmults:
                        print(f"WARNING: multiplicity for {frname} not defined! Set it to 1 by default.")
                        self.spinmults[frname] = 1
                    else:
                        print(f"Set multiplicity for {frname} to {self.spinmults[frname]}")
                else:
                    print(f"{frname} is a known fragment.")
        else:
            raise ValueError("input of charge or multiplicity not accepted")

    def computeNetCharge(self):
        netcharge = 0
        for frname in self.fragnames:
            if frname in self.newresiduenames:
                netcharge += self.charges[frname]
        self.netcharge = netcharge

    def getFragmentAtomIndex(self):
        """Get the exact atom labels in the original system for each new fragments."""
        fragresidueatomidxs = []
        for i, frag in enumerate(self.fragmols):
            frag:ob.OBMol
            fragresidueatomidx = [frag.GetAtom(j).GetId() + 1 for j in range(1, frag.NumAtoms() + 1)]
            fragresidueatomidxs.append(fragresidueatomidx)
        return fragresidueatomidxs

    def updateAtomLabels(self):
        r"""
        Change the atom label in pdb to the standard amber format. The atom label in the mainpdb may be incorrect, which may cause problems when running tleap. Only the NEW residues will be updated.
        """
        # get standard atom names for each individual residues.
        res_aname_dict = {}
        for fragpdb in self.newfragpdbs:
            atomnames = []
            f = open(fragpdb, "r")
            for line in f:
                if not (line.startswith("ATOM") or line.startswith("HETATM")):
                    continue
                atomidx = line[6:11]
                atomname = line[12:16]
                resname = line[17:20]
                atomnames.append((int(atomidx), atomname))
            f.close()
            res_aname_dict[resname] = atomnames

        # reduce complexity by pre determine the line index for each atom
        natom = self.mol_obmol.NumAtoms()
        atomlinelocations = [-1 for i in range(natom)]
        with open(self.mainfilename, "r") as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atomlinelocations[int(line[6:11])-1] = i

        # generate atom labels sequentially into a list
        fragatomidxs = self.getFragmentAtomIndex()
        for fragname, fragatomidx in zip(self.fragnames, fragatomidxs):
            if fragname not in res_aname_dict:
                continue
            fragatomlabel = res_aname_dict[fragname]
            for aid, (afid, aflb) in zip(fragatomidx, fragatomlabel):
                targetline = lines[atomlinelocations[aid-1]]
                targetline = targetline[0:12] + aflb + targetline[16:]
                lines[atomlinelocations[aid-1]] = targetline
        # Determine the last line for each fragment. 'TER\n' will be added after these lines
        fragmentlastlineindex = []
        for fragatomidx in fragatomidxs:
            fragmentlastlineindex.append(max([atomlinelocations[aid-1] for aid in fragatomidx]))
        for lineidx in fragmentlastlineindex:
            if lines[lineidx+1].find("TER") == -1:
                lines[lineidx] += "TER\n"
        with open(self.mainfilename, "w") as f:
            for i, line in enumerate(lines):
                f.write(line)

    def getResLabel(self, i:int):
        """generate the residue labels. UAA,UAB,UAC...UAZ,UBA,UBB... Can support up to 676 different labels!"""
        A2Z = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        i1 = i // 26**2
        i2 = i % 26**2 // 26
        i3 = i % 26
        return "U" + A2Z[i2] + A2Z[i3]

    def checkProteinFragment(self, fragment_obmol:ob.OBMol, fragment_index):
        for j in range(fragment_obmol.NumResidues()):
            if fragment_obmol.GetResidue(j).GetName() not in self.aminoacidresidues:
                raise NotImplementedError(f"""
                The fragment {fragment_index} contains {fragment_obmol.NumResidues()} residues, of which the {j}th residue {fragment_obmol.GetResidue(j).GetName()} is not a canonical amino acid residue. Autosolvate currently does not support automatically fitting a force field for a residue. If the original system does not contain this residue, please check whether the structure of the system is reasonable.
                """)

    def generate_fragment_key(self, frag:ob.OBMol):
        """Generate a unique key for a molecule. The atom sequence matters. Current it is accomplished by computing the smiles."""
        pybel.Molecule(ob.OBMol(frag)).write("smi", "temp.smi", overwrite = True)
        with open("temp.smi", "r") as f:
            smi = f.read()
        return smi + "-"+"".join([a.GetType() for a in ob.OBMolAtomIter(frag)])

    def getFragments(self):
        r"""
        Process a structure file contains multiple components. Generates pdb files for each individual components with their resname modified. Generates a pdb file with new residue name 
        """
        self.mol_obmol:ob.OBMol
        canonicalsmiles_resname_dict = {}
        fragments = self.mol_obmol.Separate()
        fragi = -1
        for i, fragment_obmol in enumerate(fragments):
            fragment_obmol:ob.OBMol
            # There may be a peptide in the xyz file. Peptides don't need to be fitted.
            if fragment_obmol.NumResidues() > 1:  
                self.checkProteinFragment(fragment_obmol, i)
                fragmentname = "-".join([fragment_obmol.GetResidue(j).GetName() for j in range(fragment_obmol.NumResidues())])
                self.fragmols.append(ob.OBMol(fragment_obmol))
                self.fragnames.append(fragmentname)
                continue
            fragres = fragment_obmol.GetResidue(0)
            fragresname = fragres.GetName()
            # This is a water molecule or an amino acid moleucle.
            if fragresname == "HOH" or fragresname in self.aminoacidresidues:
                self.fragmols.append(ob.OBMol(fragment_obmol))
                self.fragnames.append(fragresname)
                continue
            # This is not a regular molecule. We need to generate the forcefield for it.
            # The atom order for each same fragment should be the same. This is guaranteed by computing smiles
            smi = self.generate_fragment_key(fragment_obmol)
            # Generate the new Residue names
            # Fragments with the same canonical smile will not be assigned to different label
            newresidueflag = False
            if smi not in canonicalsmiles_resname_dict:  
                newresidueflag = True
                if fragresname.startswith("U"):
                    fragi += 1
                    canonicalsmiles_resname_dict[smi] = self.getResLabel(fragi)
                else:
                    canonicalsmiles_resname_dict[smi] = fragment_obmol.GetResidue(0).GetName()
            # Set the new residue name to the target fragment and the corresponding residue in original pdb
            resname = canonicalsmiles_resname_dict[smi]
            if fragresname.startswith("U"):
                fragres.SetName(resname)
            for j in range(1, fragment_obmol.NumAtoms()+1):
                self.mol_obmol.GetAtom(fragment_obmol.GetAtom(j).GetId()+1).GetResidue().SetName(resname)
                fragment_obmol.GetAtom(j).GetResidue().SetName(resname)
            # Output the individual fragments to enable antechamber charge fitting
            if newresidueflag:
                fragmol = pybel.Molecule(ob.OBMol(fragment_obmol))
                fragmol.write("pdb", self.basename + "-" + resname.lower() + ".pdb", overwrite = True)
                self.newfragpdbs.append(self.basename + "-" + resname.lower() + ".pdb")
                self.newfragmols.append(ob.OBMol(fragment_obmol))
                self.newresiduenames.append(resname)
            self.fragmols.append(fragment_obmol)
            self.fragnames.append(resname)

        if os.path.exists("temp.smi"):
            os.remove("temp.smi")

    def getFrcmod(self):
        cmd4 = f"$AMBERHOME/bin/parmchk2 -i {self.basename}.mol2 -f mol2 -o {self.basename}.frcmod"
        if self.srun_use:
            cmd4 = 'srun -n 1 ' + cmd4
        subprocess.run(cmd4, shell=True)

    def writeTleapcmdMultipleFragment(self):
        r"""
        Write tleap file for Multiple fragment file

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        f = open(f"{self.pathname}/leap_{self.name}.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
        for fragmentpdb in self.newfragpdbs:
            fragname, ext = os.path.splitext(fragmentpdb)
            solvent_frcmod_path = os.path.join(fragname + ".frcmod")
            solvent_prep_path   = os.path.join(fragname + ".prep")
            solvent_mol2_path   = os.path.join(fragname + ".mol2")
            solvent_lib_path    = os.path.join(fragname + ".lib")
            f.write("loadamberparams " + solvent_frcmod_path + "\n")
            for command, fpath in zip(["loadamberprep", "loadoff", "loadmol2"], [solvent_prep_path, solvent_lib_path, solvent_mol2_path]):
                if os.path.exists(fpath):
                    f.write(command + " " + fpath + "\n")
                    break
            else:
                for i in range(10):         
                    print(solvent_prep_path, solvent_mol2_path, solvent_lib_path)
                raise FileNotFoundError(f"ERROR: no amberprep or library file for fragment {fragname}")
        f.write("\n")
        f.write("SYS = loadpdb " + self.mainfilename + "\n")
        f.write("check SYS\n")
        f.write("\n")
        f.write(f"saveoff SYS {self.basename}.lib\n")
        f.write(f"saveamberparm SYS {self.basename}.prmtop {self.basename}.inpcrd\n")
        f.write(f"savepdb SYS {self.basename}.pdb\n")
        f.write(f"savemol2 SYS {self.basename}.mol2 1\n")
        f.write("quit\n")
        f.close()

    def buildAmberParamsForAll(self):
        self.reorderPDB(self.mainfilename)
        self.getFragments()
        self.formatPDB(self.mainfilename)
        for fragpdb in self.newfragpdbs:
            self.formatPDB(fragpdb)
        self.checkParams()
        self.computeNetCharge()
        pybel.Molecule(ob.OBMol(self.mol_obmol)).write("pdb", self.mainfilename, overwrite=True)
        
        for i, (resname, pdbname) in enumerate(zip(self.newresiduenames, self.newfragpdbs)):
            if self.charges[resname] == 0 and self.spinmults[resname] == 1 and self.pre_optimize_fragments: # pre-optimize each fragment before proceeding to antechamber
                bname = os.path.splitext(pdbname)[0]
                os.rename(pdbname, bname+"1.pdb")
                subprocess.run(f"obminimize -ff gaff {bname}1.pdb > {pdbname}", shell = True)
                os.remove(bname+"1.pdb")
            mol = pybel.readfile("pdb", pdbname).__next__()
            mname = os.path.splitext(os.path.basename(pdbname))[0]
            inst = AmberParamsBuilder(
                pdbname,
                name = mname,
                resname = mol.OBMol.GetResidue(0).GetName(),
                charge = self.charges[resname],
                spinmult = self.spinmults[resname],
                charge_method = self.charge_method,
            )
            inst.build()
        self.updateAtomLabels()
        self.writeTleapcmdMultipleFragment()
        cmd =f"tleap -s -f {self.pathname}/leap_{self.name}.cmd > {self.pathname}/leap_{self.name}.log"
        subprocess.call(cmd, shell=True)
        self.checkLeapOutput(f"{self.pathname}/leap_{self.name}.log")
        # self.getFrcmod()
        if self.deletefiles:
            self.clean_cwd()

    def checkLeapOutput(self, logfile):
        with open(logfile, "r") as f:
            for line in f:
                if line.find("FATAL") != -1:
                    raise AssertionError("Error in running tleap. " + line)

    def clean_cwd(self):
        cwd = os.getcwd()
        basename = os.path.splitext(os.path.basename(self.mainfilename))[0]
        pathname = os.path.dirname(self.mainfilename)
        if os.path.exists(pathname):
            os.chdir(pathname)
        fnames = os.listdir(os.getcwd())
        for fname in fnames:
            if fname.startswith("ANTECHAMBER"):
                os.remove(fname)
        if os.path.exists("ATOMTYPE.INF"):
            os.remove("ATOMTYPE.INF")
        for ufres in self.newresiduenames:
            print(ufres, os.listdir(os.getcwd()))
            ufres = ufres.lower()
            if os.path.exists(f"leap_{basename}-{ufres}_savelib.log"):
                os.remove(f"leap_{basename}-{ufres}_savelib.log")
            if os.path.exists(f"leap_{basename}-{ufres}.cmd"):
                os.remove(f"leap_{basename}-{ufres}.cmd")
            if os.path.exists(f"leap_{basename}-{ufres}.log"):
                os.remove(f"leap_{basename}-{ufres}.log")
            for suffix in [".frcmod", ".lib", ".mol2", "_nonprot.pdb", "_renum.txt", "_sslink", "-format.log", ".pdb", ".inpcrd", ".prmtop"]:
                if os.path.exists(f"{basename}-{ufres}{suffix}"):
                    os.remove(f"{basename}-{ufres}{suffix}")
        for suffix in [".frcmod", ".lib", ".mol2", "_nonprot.pdb", "_renum.txt", "_sslink", "-format.log"]:
            if os.path.exists(f"{basename}{suffix}"):
                os.remove(f"{basename}{suffix}")
        if os.path.exists(f"leap_{basename}.log"):
            os.remove(f"leap_{basename}.log")
        if os.path.exists(f"leap_{basename}.cmd"):
            os.remove(f"leap_{basename}.cmd")
        if os.path.exists(f"leap.log"):
            os.remove(f"leap.log")
        for f in ["sqm.in","sqm.out","sqm.pdb","stdout_nonprot.pdb","stdout_renum.pdb","stdout_sslink"]:
            if os.path.exists(f):
                os.remove(f)
        os.chdir(cwd)
        
class MulticomponentSolventBoxBuilder(solventBoxBuilder):
    def __init__(self, xyzfile: str, 
    slu_charge=0, slu_spinmult=1, slu_count=1, 
    solvent="water", solvent_frcmod="", solvent_off="",
    slv_xyz="", slv_generate=False, slv_count=210*8, cube_size = 54, closeness = 0.8,
    charge_method="resp", outputFile="", pre_optimize_fragments = False, deletefiles=False,
    srun_use=False, gaussianexe=None, gaussiandir=None, amberhome=None, use_terachem = False): 

        super().__init__(xyzfile,
        solvent = solvent, slu_netcharge=slu_charge, slu_spinmult=slu_spinmult, 
        cube_size = cube_size, charge_method = charge_method, outputFile = outputFile, srun_use=srun_use,
        gaussianexe = gaussianexe, gaussiandir=gaussiandir, amberhome=amberhome, use_terachem=use_terachem,
        closeness=closeness, solvent_off = solvent_off, solvent_frcmod=solvent_frcmod, slu_count=slu_count, slv_count=slv_count, slv_generate=slv_generate, slv_xyz=slv_xyz)

        self.pre_optimize_fragments = pre_optimize_fragments
        self.deletefiles = deletefiles
        self.slu_charge = slu_charge
        self.newresiduenames = []
        self.newfragpdbs = []

    def inputCheck(self):
        if self.amberhome == None:
            print("WARNING: Amber home directory is not specified in input options")
            print("WARNING: Checking AMBERHOME environment variable...")
            cmd = ["echo", "$AMBERHOME"]
            print(cmd)
            proc=subprocess.Popen(cmd, universal_newlines=True, stdout=subprocess.PIPE)
            output=proc.communicate()[0]
            if len(output)<2:
                raise ModuleNotFoundError("ERROR: AMBERHOME not defined! Please set up AMBERHOME environment or specify through command line option -a")
            else:
                print("AMBERHOME detected: ", output)
                self.amberhome = output
        else:
            print("AMBERHOME path provided from input: ", self.amberhome)
            print("Validating path...")
            if not os.path.isdir(self.amberhome):
                raise FileNotFoundError("ERROR: provided AMBERHOME path does not exist! Exiting...")
            print("Exporting AMBERHOME environment variable:")
            cmd = "export AMBERHOME=" + self.amberhome
            print(cmd)
            subprocess.call(cmd, shell=True)
            print("AMBERHOME environment variable export finished.")
        # check whether requested custom solvent is available
        if self.solvent not in amber_solv_dict.keys() and self.solvent not in custom_solv_dict.keys():
            print("Requested solvent name not contained in AutoSolvate.")
            print("Checking available of custom solvent .frcmod and .off files")
            if not self.slv_generate:
                if len(self.solvent_frcmod) == 0:
                    raise FileNotFoundError("ERROR!Custom solvent .frcmod file is not provided!")
                elif len(self.solvent_off) == 0:
                    raise FileNotFoundError("ERROR!Custom solvent .off library file is not provided!")
                elif not os.path.exists(self.solvent_frcmod):
                    raise FileNotFoundError("ERROR!Custom solvent .frcmod file ", self.solvent_frcmod, "does not exist!")
                elif not os.path.exists(self.solvent_off):
                    raise FileNotFoundError("ERROR!Custom solvent .off library file ", self.solvent_frcmod, "does not exist!")
                else:
                    self.is_custom_solvent = True
            elif self.slv_generate:
                if len(self.solvent_frcmod) == 0:
                    print("Custom solvent .frcmod file is not provided! Will generate from GAFF instead.")
                elif len(self.solvent_off) == 0:
                    print("Custom solvent .off library file is not provided! The box will generate from packmol instead.")
                elif not os.path.exists(self.solvent_frcmod):
                    print("Custom solvent .frcmod file ", self.solvent_frcmod, "does not exist! Will generate from GAFF instead.")
                elif not os.path.exists(self.solvent_off):
                    print("Custom solvent .off library file ", self.solvent_frcmod, "does not exist! The box will generate from packmol instead.")
                else:
                    print("Solvent parameter generation is disabled as .frcmod and .off are provided.")
                    self.is_custom_solvent = True
                    self.slv_generate = False
            if self.slv_generate:
                if not self.slv_xyz or not os.path.exists(self.slv_xyz):
                    raise FileNotFoundError("ERROR!The coordinate file for custom solvent must be provided when generate parameter from forcefield.")
                cmd = "which packmol"
                proc=subprocess.Popen(cmd, shell = True, universal_newlines=True, stdout=subprocess.PIPE)
                output=proc.communicate()[0]
                if not output:
                    raise ModuleNotFoundError("ERROR! Cannot find packmol in $PATH. packmol must be installed when using custom solvent and multiple solute.")
            
        else:
            print("Solvent parameter generation is disabled as requested solvent is contained in Autosolvate.")
            self.slv_generate = False
        

        if len(self.solute.atoms) > 100:
            print("The solute molecule has over 100 atoms. Will adjust -pl to 15")
   
    def writeTleapcmd_custom_solvated(self):
        r"""
        Write tleap file for custom solvents like CH3CN

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        if self.slv_generate:
            solvPrefix = self.solvent
            solvent_frcmod = {self.slv_name}+"frcmod"
            solvent_frcmod_path = os.path.join(os.getcwd(), solvent_frcmod)
            solvent_prep = {self.slv_name}+"prep"
            solvent_prep_path = os.path.join(os.getcwd(), solvent_prep)
            solvent_mol2 = {self.slv_name}+"mol2"
            solvent_mol2_path = os.path.join(os.getcwd(), solvent_mol2)
            solvent_lib = {self.slv_name}+"lib"
            solvent_lib_path = os.path.join(os.getcwd(), solvent_lib)
        else:
            solvPrefix = custom_solv_dict[self.solvent]
            solvent_frcmod = solvPrefix+'.frcmod'
            solvent_frcmod_path = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data',solvPrefix,solvent_frcmod))

            solvent_prep = solvPrefix+'.prep'
            solvent_prep_path = pkg_resources.resource_filename('autosolvate', 
                    os.path.join('data',solvPrefix,solvent_prep))
            solvent_lib_path = ""
            solvent_mol2_path = ""

        f = open("leap_packmol_solvated.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n") # This will load the Ions. Neccessary
        f.write("loadamberparams " + solvent_frcmod_path + "\n")
        for command, fpath in zip(["loadamberprep", "loadoff", "loadmol2"], [solvent_prep_path, solvent_lib_path, solvent_mol2_path]):
            if not fpath:
                continue
            if os.path.exists(fpath):
                f.write(command + " " + fpath + "\n")
                break
        else:
            print("ERROR: no solvent amberprep or library file!")
            exit()
        for pdbname in self.newfragpdbs:
            bname = os.path.splitext(os.path.basename(pdbname))[0]
            f.write(f"loadamberparams {bname}.frcmod\n")
            f.write(f"loadoff {bname}.lib\n")
        f.write("\n")
        f.write("SYS = loadpdb " + solvPrefix + "_solvated.processed.pdb\n")
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
        f.write("saveamberparm SYS "+str(self.outputFile)+".prmtop "+str(self.outputFile)+".inpcrd #Save AMBER topology and coordinate files\n")
        f.write("savepdb SYS "+str(self.outputFile)+".pdb\n")
        f.write("quit\n")
        f.close()

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
        if self.solvent in amber_solv_dict:
            print("Now add pre-equlibrated solvent box to the solute")
            f = open("leap_add_solventbox.cmd","w")
            f.write("source leaprc.protein.ff14SB\n")
            f.write("source leaprc.gaff\n")
            f.write("source leaprc.water.tip3p\n")
            f.write(str(amber_solv_dict[str(self.solvent)][0]))
            for pdbname in self.newfragpdbs:
                bname = os.path.splitext(os.path.basename(pdbname))[0]
                f.write(f"loadamberparams {bname}.frcmod\n")
                f.write(f"loadoff {bname}.lib\n")
            f.write(f"mol=loadmol2 {self.slu_name}.mol2\n")
            f.write("check mol\n")
            f.write("solvatebox mol " + (str(amber_solv_dict[str(self.solvent)][1])) + str(self.slu_pos) + " iso "+str(self.closeness)+"  #Solvate the complex with a cubic solvent box\n") 
            # Notice that we want to add the ion after solvation because we don't want the counter ion to be too close to solute
            if self.slu_netcharge != 0:
                if self.slu_netcharge > 0:
                    ion = 'Cl-'
                else:
                    ion = 'Na+'
                f.write("addIons2 mol " + ion + " 0\n")
                f.write("check mol\n")
            f.write("check mol\n")
            f.write("savepdb mol " + str(self.outputFile) + ".pdb\n")
            f.write("saveamberparm mol " + str(self.outputFile) + ".prmtop " + str(self.outputFile) + ".inpcrd\n")
            f.write("quit\n")
            f.close()

    def writeTleapcmd_add_solvent_custom(self):
        r"""
        Write tleap input file to add solvent from user provided .off and .frcmod files

        Parameters
        ----------
        None

        Returns
        -------
        None
        """
        print("Now add custom pre-equlibrated solvent box to the solute")
        f = open("leap_add_solventbox.cmd","w")
        f.write("source leaprc.protein.ff14SB\n")
        f.write("source leaprc.gaff\n")
        f.write("source leaprc.water.tip3p\n")
        f.write("loadoff " + self.solvent_off + "\n")
        f.write("loadamberparams " + self.solvent_frcmod + "\n")
        for pdbname in self.newfragpdbs:
            bname = os.path.splitext(os.path.basename(pdbname))[0]
            f.write(f"loadamberparams {bname}.frcmod\n")
            f.write(f"loadoff {bname}.lib\n")
        f.write(f"mol=loadmol2 {self.slu_name}.mol2\n")
        f.write("check mol\n")
        f.write("solvatebox mol " + self.solvent + " " 
                + str(self.slu_pos) + " iso 0.8  #Solvate the complex with a cubic solvent box\n") 
        # Notice that we want to add the ion after solvation because we don't want the counter ion to be too close to solute
        if self.slu_netcharge != 0:
            if self.slu_netcharge > 0:
                ion = 'Cl-'
            else:
                ion = 'Na+'
            f.write("addIons2 mol " + ion + " 0\n")
            f.write("check mol\n")
        f.write("check mol\n")
        f.write("savepdb mol " + str(self.outputFile) + ".pdb\n")
        f.write("saveamberparm mol " + str(self.outputFile) + ".prmtop " + str(self.outputFile) + ".inpcrd\n")
        f.write("quit\n")
        f.close()

    def build(self):
        if check_multicomponent(self.xyz):
            self.solutebuilder = MulticomponentParamsBuilder(
                xyzfile=self.xyz,
                name = "",
                resname="",
                charge = self.slu_netcharge,
                spinmult=self.slu_spinmult,
                charge_method=self.charge_method,
                outputFile=self.outputFile,
                pre_optimize_fragments=self.pre_optimize_fragments,
                srun_use=self.srun_use,
                gaussianexe=self.gaussian_exe,
                gaussiandir=self.gaussian_dir,
                amberhome = self.amberhome,
                deletefiles=False
            )
            self.solutebuilder.buildAmberParamsForAll()
            self.slu_netcharge = self.solutebuilder.netcharge
            self.newresiduenames = self.solutebuilder.newresiduenames
            self.newfragpdbs = self.solutebuilder.newfragpdbs   
        else:
            self.solutebuilder = AmberParamsBuilder(
                xyzfile= self.xyz,
                name = self.slu_name,
                resname = "SLU",
                charge = self.slu_netcharge,
                spinmult= self.slu_spinmult,
                charge_method=self.charge_method,
                outputFile="solutegen.out",
                srun_use=self.srun_use,
                gaussianexe=self.gaussian_exe,
                gaussiandir=self.gaussian_dir,
                amberhome=self.amberhome
            )
            self.solutebuilder.build()
            self.newresiduenames = [self.solutebuilder.resname, ]
            self.newfragpdbs = [self.solutebuilder.xyz, ]

        if self.slv_generate:
            solventbuilder = AmberParamsBuilder(
                xyzfile = self.slv_xyz,
                name = self.solvent,
                resname = "SLV",
                charge = 0,
                spinmult = 1,
                charge_method=self.charge_method,
                outputFile="solventgen.out",
                srun_use=self.srun_use,
                gaussianexe=self.gaussian_exe,
                gaussiandir=self.gaussian_dir, 
                amberhome=self.amberhome
            )
            solventbuilder.build()
            
        self.createAmberParm()
        print("The script has finished successfully")





if __name__ == "__main__":
    inst = MulticomponentParamsBuilder("PAHs.pdb", deletefiles=True)
    inst.buildAmberParamsForAll()
