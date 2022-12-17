#--------------------------------------------------------------------------------------------------#
# multicomponent.py. 
# Description: 
#   This module can handle structure files containing multiple molecules.
#   1. Generates the lib and frcmod file for each separate molecule when it is not h2o or amino acids.
#   2. Generates the lib, prmtop and inpcrd files of this whole structure.
# author: Fangning Ren (2022-12-17) 
# path: autosolvate/multicomponent.py
#--------------------------------------------------------------------------------------------------#
import getopt, sys, os
import subprocess
from openbabel import pybel
from openbabel import openbabel as ob

from autosolvate import AmberParamsBuilder


class MulticomponentParamsBuilder():
    def __init__(self, xyzfile: str, 
    name="", resname="", charge=0, spinmult=1, 
    charge_method="resp", outputFile="", srun_use=False, gaussianexe=None, gaussiandir=None, amberhome=None): 
        """
        The charge and spinmultiplicity for each fragment should be sequentially defined by the user. 
        The length of these array should equal to the number of fragments.
        Otherwise they will be set to zero and one by default. 
        """
        self.names = name
        self.resnames = resname
        self.charges = charge
        self.spinmults = spinmult

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
            f = open(self.basename + "-origin" + ext, "w")
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

        self.aminoacidresidues = set(["ALA", "ARG", "ASN", "ASP", "CYS", "GLU", "GLN", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"])

        self.getPDB()
        self.mol_obmol = pybel.readfile("pdb", self.mainfilename).__next__().OBMol
        self.fragnames = []
        self.fragmols = []
        self.newfragmols = []
        self.newfragpdbs = []
        self.newresiduenames = []

    def getPDB(self):
        basename, ext = os.path.splitext(self.mainfilename)
        ext = ext[1:]
        if ext != "pdb":
            subprocess.run(f"obabel -i {ext} {self.mainfilename} -o pdb > {basename}.pdb", shell = True)
        self.mainfilename = basename + ".pdb"

    def formatPDB(self, pdbfile:str):
        mainname, ext = os.path.splitext(pdbfile)
        os.rename(pdbfile, mainname + "-original.pdb")
        subprocess.run(f"pdb4amber {mainname}-original.pdb > {pdbfile}", shell = True)
        os.remove(f"{mainname}-original.pdb")

    def checkParams(self):
        n_frag = len(self.fragmols)
        if isinstance(self.charges, int) or isinstance(self.spinmults, int):
            print("WARNING: All charge and multiplicity are settled to 0 and 1")
            self.charges = [0 for i in range(n_frag)]
            self.spinmults = [1 for i in range(n_frag)]
            return
        if len(self.charges) != len(self.fragmols) or len(self.spinmults) != len(self.fragmols):
            raise IndexError(f"The length of the provided charge or multiplicity array are different from the number of segments in the system.")

    def getFragmentAtomIndex(self):
        r"""
        Get the exact atom labels in the original system for each new fragments.

        Parameters
        ---------
        mainpdb:

        The input pdb with multiple fragments

        Returns
        ---------
        list[list[]]

        """
        fragresidueatomidxs = []
        for i, frag in enumerate(self.fragmols):
            frag:ob.OBMol
            fragresidueatomidx = [frag.GetAtom(j).GetId() + 1 for j in range(1, frag.NumAtoms() + 1)]
            fragresidueatomidxs.append(fragresidueatomidx)
        return fragresidueatomidxs

    def updateAtomLabels(self):
        r"""
        Change the atom label in pdb to the standard amber format. The atom label in the mainpdb may be incorrect, which may cause problems when running tleap. Only the NEW residues will be updated.
        Parameters
        ---------
        mainpdb:

        The input pdb with multiple fragments

        fragpdbs:

        a list of tleap generated pdb files that represents different fragment of the system
        
        Returns
        ---------
        None

        The atom label in mainpdb will be relabeled

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

    def getFragments(self):
        r"""
        Process a structure file contains multiple components. Generates pdb files for each individual components with their resname modified. Generates a pdb file with new residue name 
            
        Parameters
        ---------
        filename:

        The input file with multiple fragments
        
        Returns
        ---------
        List[str]

        List of the pdb file names of each individual components

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
            pybel.Molecule(ob.OBMol(fragment_obmol)).write("smi", "temp.smi", overwrite = True)
            with open("temp.smi", "r") as f:
                smi = f.read()
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
        f.write("quit\n")
        f.close()

    def buildAmberParamsForAll(self):
        self.getFragments()
        self.formatPDB(self.mainfilename)
        for fragpdb in self.newfragpdbs:
            self.formatPDB(fragpdb)
        self.checkParams()
        pybel.Molecule(ob.OBMol(self.mol_obmol)).write("pdb", self.mainfilename, overwrite=True)
        
        for i, pdbname in enumerate(self.newfragpdbs):
            mol = pybel.readfile("pdb", pdbname).__next__()
            mname = os.path.splitext(os.path.basename(pdbname))[0]
            inst = AmberParamsBuilder(
                pdbname,
                name = mname,
                resname = mol.OBMol.GetResidue(0).GetName(),
                charge = self.charges[i],
                spinmult = self.spinmults[i],
                charge_method = "bcc",
            )
            inst.build()
        self.updateAtomLabels()
        self.writeTleapcmdMultipleFragment()
        cmd =f"tleap -s -f {self.pathname}/leap_{self.name}.cmd > {self.pathname}/leap_{self.name}.log"
        subprocess.call(cmd, shell=True)
        self.checkLeapOutput(f"{self.pathname}/leap_{self.name}.log")

    def checkLeapOutput(self, logfile):
        with open(logfile, "r") as f:
            for line in f:
                if line.find("FATAL") != -1:
                    raise AssertionError("Error in running tleap. " + line)


if __name__ == "__main__":
    inst = MulticomponentParamsBuilder("PAHs.pdb")
    inst.buildAmberParamsForAll()