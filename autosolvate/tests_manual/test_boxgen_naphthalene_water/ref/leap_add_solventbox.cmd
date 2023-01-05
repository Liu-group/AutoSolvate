source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p
 loadamberparams solute.frcmod
mol=loadmol2 solute.mol2
check mol
solvatebox mol TIP3PBOX 27.0 iso 0.8  #Solvate the complex with a cubic water box
check mol
savepdb mol neutral.pdb
saveamberparm mol neutral.prmtop neutral.inpcrd
quit
