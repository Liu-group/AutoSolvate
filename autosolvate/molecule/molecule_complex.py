import getopt
import sys
import os
import shutil
import subprocess
import pkg_resources
import glob
import logging

from typing import Any, List, Tuple
from dataclasses import dataclass, field, asdict

from openbabel import pybel
from openbabel import openbabel as ob

from ..Common import *
from ..utils import *
from .molecule import *

from logging import DEBUG, INFO, WARN, CRITICAL
logging.basicConfig(level = INFO, force = True, handlers=[])

logger              = logging.getLogger(name = "MoleculeComplex")
# output_handler      = logging.FileHandler(filename = "log.txt", mode = "a", encoding="utf-8")
output_handler      = logging.StreamHandler()
output_formater     = logging.Formatter(fmt = '%(asctime)s %(name)s %(levelname)s: %(message)s', datefmt="%H:%M:%S")
output_handler.setFormatter(output_formater)
if len(logger.handlers) == 0:
    logger.addHandler(output_handler)


class MoleculeComplex(System):

    _SUPPORT_INPUT_FORMATS = ['pdb', 'xyz']

    def __init__(self, 
                xyzfile:        str = "",
                charges:        int = 0,
                multiplicities: int = 1,
                name:           str = "", 
                residue_name:   str = "SYS",
                folder:         str = ".",
                reorder_pdb:    bool = False,
                ):
        self.name           = process_system_name(name, xyzfile, support_input_format=MoleculeComplex._SUPPORT_INPUT_FORMATS)   
        self.folder         = os.path.abspath(folder)
        self.charges        = charges
        self.multiplicity   = multiplicities
        self.spinmults      = multiplicities
        self.residue_name   = "SYS" if not residue_name else residue_name
        self.number         = 0
        self.read_coordinate(xyzfile)
        super(MoleculeComplex, self).__init__(name = self.name)

        self.mol_obmol          = pybel.readfile("pdb", self.pdb).__next__().OBMol
        self.fragresiduenames   = []    # name of each fragment
        self.fragmols           = []    # ob.OBMol object for each fragment
        self.newfragmols        = []    # ob.OBMol object for new fragments that are not known residues, ligands, or water
        self.newfragpdbs        = []    # pdb file for these new fragments
        self.newresiduenames    = []    # residue name for these new fragments

        self.newmolecules       = []    # list of Molecule objects
        self.netcharge          = 0     # net charge of the system
        self.multiplicity       = 1     # multiplicity of the system

        if reorder_pdb:
            reorderPDB(self.pdb, self.pdb)
        self.build_molecules()

    def read_coordinate(self, fname:str):
        ext = os.path.splitext(fname)[-1][1:]
        setattr(self, ext, fname)
        for e in Molecule._SUPPORT_INPUT_FORMATS:
            if e == ext:
                continue
            nname = os.path.splitext(fname)[0] + "." + e
            subprocess.run(f"obabel -i {ext} {fname} -o {e} > {nname}", shell = True)
            setattr(self, e, nname)

    def check_protein_fragment(self, fragment_obmol:ob.OBMol, fragment_index = 0):
        for j in range(fragment_obmol.NumResidues()):
            if fragment_obmol.GetResidue(j).GetName() not in self.aminoacidresidues:
                logger.critical(f"The fragment {fragment_index} contains {fragment_obmol.NumResidues()} residues")
                logger.critical(f"of which the {j}th residue {fragment_obmol.GetResidue(j).GetName()} is not a canonical amino acid residue. ")
                logger.critical(f"Autosolvate currently does not support automatically fitting a force field for a residue. ")
                logger.critical(f"If the original system does not contain this residue, please check whether the structure of the system is reasonable.""")
                raise NotImplementedError("see the log for details")

    def generate_fragment_key(self, frag:ob.OBMol, atomorder = True, smiles = True):
        """Generate a unique key for a molecule. The atom sequence matters. Current it is accomplished by computing the smiles."""
        aord = "".join([a.GetType() for a in ob.OBMolAtomIter(frag)]) if atomorder else ""
        if smiles:
            pybel.Molecule(ob.OBMol(frag)).write("smi", "temp.smi", overwrite = True)
            with open("temp.smi", "r") as f:
                smi = f.read()
        else:
            smi = ""
        return smi + "-" + aord

    def get_residue_label(self, i:int):
        """generate the residue labels. UAA,UAB,UAC...UAZ,UBA,UBB... Can support up to 676 different labels!"""
        A2Z = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
        i1 = i // 26**2
        i2 = i % 26**2 // 26
        i3 = i % 26
        return "U" + A2Z[i2] + A2Z[i3]

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
                self.check_protein_fragment(fragment_obmol, i)
                fragmentname = "-".join([fragment_obmol.GetResidue(j).GetName() for j in range(fragment_obmol.NumResidues())])
                self.fragmols.append(ob.OBMol(fragment_obmol))
                self.fragresiduenames.append(fragmentname)
                logger.info(f"Fragment {i} is a peptide with {fragment_obmol.NumResidues()} residues. ")
                continue
            fragres = fragment_obmol.GetResidue(0)
            fragresname = fragres.GetName()
            # This is a water molecule or an amino acid moleucle.
            if fragresname == "HOH" or fragresname in AMINO_ACID_RESIDUES:
                self.fragmols.append(ob.OBMol(fragment_obmol))
                self.fragresiduenames.append(fragresname)
                logger.info(f"Fragment {i} is a known molecule with res name {fragresname}. ")
                continue
            # This is not a regular molecule. We need to generate the forcefield for it.
            # The atom order for each same fragment should be the same. This is guaranteed by computing smiles
            smi = self.generate_fragment_key(fragment_obmol)
            # Generate the new Residue names
            # Fragments with the same canonical smile will not be assigned to different label
            newresidueflag = False
            if smi not in canonicalsmiles_resname_dict:  
                newresidueflag = True
                if fragresname.startswith("U"): # This new fragment has the default residue name assigned by openbabel
                    fragi += 1
                    canonicalsmiles_resname_dict[smi] = self.get_residue_label(fragi)
                else:   # This new fragment has a unique residue name assigned by the user
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
                fragmol.write("pdb", self.reference_name + "-" + resname.lower() + ".pdb", overwrite = True)
                self.newfragpdbs.append(self.reference_name + "-" + resname.lower() + ".pdb")
                self.newfragmols.append(ob.OBMol(fragment_obmol))
                self.newresiduenames.append(resname)
                logger.info(f"Fragment {i} is a new molecule with res name {resname}. Update the term list.")
            else:
                logger.info(f"Fragment {i} is a known molecule with res name {resname}. ")
            self.fragmols.append(fragment_obmol)
            self.fragresiduenames.append(resname)

        if os.path.exists("temp.smi"):
            os.remove("temp.smi")

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
        with open(self.pdb, "r") as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if not (line.startswith("ATOM") or line.startswith("HETATM")):
                continue
            atomlinelocations[int(line[6:11])-1] = i

        # generate atom labels sequentially into a list
        fragatomidxs = self.getFragmentAtomIndex()
        for fragname, fragatomidx in zip(self.fragresiduenames, fragatomidxs):
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
        with open(self.pdb, "w") as f:
            for i, line in enumerate(lines):
                f.write(line)


    def checkParams(self):
        """check parameters, especially for the length of charge and multiplicity"""
        logger.info(f"All fragments: {' '.join(self.fragresiduenames)}")
        logger.info(f"New fragments: {' '.join(self.newresiduenames)}")
        if isinstance(self.charges, int):
            logger.warn("All charges are set to 0")
            self.charges = {frname:0 for frname in self.fragresiduenames}
        if isinstance(self.spinmults, int):
            logger.warn("All multiplicities are set to 1")
            self.spinmults = {frname:1 for frname in self.fragresiduenames}
        if isinstance(self.charges, list):
            if len(self.charges) != len(self.fragresiduenames):
                logger.critical(f"Only charge for {len(self.charges)} fragments are provided. This file has {len(self.fragresiduenames)} fragments!")
                raise ValueError(f"Only charge for {len(self.charges)} fragments are provided. This file has {len(self.fragresiduenames)} fragments!")
            self.charges = {frname:c for (frname, c) in zip(self.fragresiduenames, self.charges)}
        if isinstance(self.spinmults, list):
            if len(self.spinmults) != len(self.fragresiduenames):
                logger.critical(f"Only multiplicity for {len(self.spinmults)} fragments are provided. This file has {len(self.fragresiduenames)} fragments!")
                raise ValueError(f"Only multiplicity for {len(self.spinmults)} fragments are provided. This file has {len(self.fragresiduenames)} fragments!")
            self.spinmults = {frname:c for (frname, c) in zip(self.fragresiduenames, self.spinmults)}
        if isinstance(self.charges, dict) and isinstance(self.spinmults, dict):
            for frname in self.fragresiduenames:
                if len(frname) > 5:
                    logger.info(f"{frname} is a peptide with standard residues.")
                elif frname in self.newresiduenames:
                    if frname[0] == "U":
                        logger.info(f"{frname} corresponds to a new fragment with auto-generated name.")
                    else:
                        logger.info(f"{frname} corresponds to a new fragment.")
                    if frname not in self.charges:
                        logger.warn(f"charge for {frname} not defined! Set it to 0 by default.")
                        self.charges[frname] = 0
                    else:
                        logger.info(f"Set charge for {frname} to {self.charges[frname]}")
                    if frname not in self.spinmults:
                        logger.warn(f"multiplicity for {frname} not defined! Set it to 1 by default.")
                        self.spinmults[frname] = 1
                    else:
                        logger.info(f"Set multiplicity for {frname} to {self.spinmults[frname]}")
                else:
                    logger.info(f"{frname} is a known fragment.")
        else:
            logger.critical("input of charge or multiplicity not accepted")
            logger.critical(f"charge: {self.charges}")
            logger.critical(f"multiplicity: {self.spinmults}")
            raise ValueError("input of charge or multiplicity not accepted")
        
    def computeNetCharge(self):
        netcharge = 0
        for frname in self.fragresiduenames:
            if frname in self.newresiduenames:
                netcharge += self.charges[frname]
        self.netcharge = netcharge
        logger.info(f"Net charge of the molecule is {netcharge}")

    def computeMultiplicity(self):
        not_paired_electrons = 0
        for frname in self.fragresiduenames:
            if frname in self.newresiduenames:
                not_paired_electrons += self.spinmults[frname] - 1
        self.multiplicity = not_paired_electrons + 1
        logger.info(f"Total multiplicity of the molecule is {self.multiplicity}")

    def build_molecules(self):
        self.getFragments()
        formatPDB(self.pdb)
        for fragpdb in self.newfragpdbs:
            formatPDB(fragpdb)
        self.updateAtomLabels()
        self.checkParams()
        self.computeNetCharge()
        self.computeMultiplicity()



        for newresiduename, newpdbname, newfragmol in zip(self.newresiduenames, self.newfragpdbs, self.newfragmols):
            logger.info(f"Create Molecule object for fragment {newresiduename}")
            charge = self.charges[newresiduename]
            spinmult = self.spinmults[newresiduename]
            mol = Molecule(newpdbname, charge, spinmult, residue_name=newresiduename, folder = self.folder)
            self.newmolecules.append(mol)
