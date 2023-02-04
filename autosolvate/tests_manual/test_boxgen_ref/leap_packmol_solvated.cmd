source leaprc.protein.ff14SB
source leaprc.gaff
source leaprc.water.tip3p
loadamberparams /opt/anaconda3/lib/python3.8/site-packages/autosolvate/data/ch3cn/ch3cn.frcmod
loadAmberPrep /opt/anaconda3/lib/python3.8/site-packages/autosolvate/data/ch3cn/ch3cn.prep
loadamberparams solute.frcmod
loadoff solute.lib

SYS = loadpdb ch3cn_solvated.processed.pdb
check SYS

# set the dimension of the periodic box
set SYS box {56, 56, 56}

saveamberparm SYS nap_neutral_MeCN.prmtop nap_neutral_MeCN.inpcrd #Save AMBER topology and coordinate files
savepdb SYS nap_neutral_MeCN.pdb
quit
