source leaprc.protein.ff14SB
source leaprc.gaff
loadamberparams solute.frcmod
SLU = loadmol2 solute.mol2
check SLU
set SLU head SLU.1.C
set SLU tail SLU.1.C9
saveoff SLU solute.lib
savepdb SLU solute.pdb
quit
