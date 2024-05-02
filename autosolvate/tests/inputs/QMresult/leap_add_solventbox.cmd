source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p
addAtomTypes {
        { "M1"  "Fe" "sp3" }
        { "Y7"  "N" "sp3" }
        { "Y8"  "N" "sp3" }
        { "Y9"  "N" "sp3" }
        { "Z1"  "N" "sp3" }
        { "Z2"  "N" "sp3" }
        { "Z3"  "N" "sp3" }
}
FE1 = loadmol2 FE1.mol2
L01 = loadmol2 L01.mol2
L11 = loadmol2 L11.mol2
L21 = loadmol2 L21.mol2
loadamberparams LG2.frcmod
loadamberparams LG0.frcmod
loadamberparams LG1.frcmod
loadamberparams frcmod.ionslm_126_opc
loadamberparams /home/sunnyxun/anaconda3/envs/autosolvate/lib/python3.8/site-packages/autosolvate-0.1.4+62.g502bde4-py3.8.egg/autosolvate/data/ch3cn/ch3cn.frcmod
loadAmberPrep /home/sunnyxun/anaconda3/envs/autosolvate/lib/python3.8/site-packages/autosolvate-0.1.4+62.g502bde4-py3.8.egg/autosolvate/data/ch3cn/ch3cn.prep
loadamberparams Fe_plus2_mcpbpy.frcmod
mol = loadpdb Fe_plus2_processed.pdb 
bond mol.1.FE mol.2.N
bond mol.1.FE mol.2.N1
bond mol.1.FE mol.3.N
bond mol.1.FE mol.3.N1
bond mol.1.FE mol.4.N
bond mol.1.FE mol.4.N1

check mol

addIons2 mol Cl- 2 
check mol
# set the dimension of the periodic box
set mol box {56.0, 56.0, 56.0}

saveamberparm mol Fe_plus2_solvated.prmtop Fe_plus2_solvated.inpcrd #Save AMBER topology and coordinate files
savepdb mol Fe_plus2_solvated.pdb
quit
