source                leaprc.protein.ff14SB       
source                leaprc.gaff                
source                leaprc.water.tip3p         
CO1 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/CO1.mol2
L01 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/L01.mol2
L11 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/L11.mol2
L21 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/L21.mol2
loadamberparams /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc_single/Co_plus3_bpy3.frcmod
TMC   = loadpdb       /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc_single/Co_plus3_dry.pdb          
bond TMC.1.CO TMC.2.N
bond TMC.1.CO TMC.2.N1
bond TMC.1.CO TMC.3.N
bond TMC.1.CO TMC.3.N1
bond TMC.1.CO TMC.4.N
bond TMC.1.CO TMC.4.N1
savepdb TMC /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc_single/Co_plus3_bpy3.pdb        
saveamberparm TMC /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc_single/Co_plus3_bpy3.prmtop /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc_single/Co_plus3_bpy3.inpcrd     
quit              
