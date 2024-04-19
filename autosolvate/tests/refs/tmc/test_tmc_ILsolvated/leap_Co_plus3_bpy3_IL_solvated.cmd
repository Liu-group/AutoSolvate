source                leaprc.protein.ff14SB       
source                leaprc.gaff                
source                leaprc.water.tip3p         
CO1 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/CO1.mol2
L01 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/L01.mol2
L11 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/L11.mol2
L21 = loadmol2 /home/fren5/AutoSolvae-update/FFgen/L21.mol2
loadamberparams /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/Co_plus3_bpy3.frcmod
TMC   = loadpdb       /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/Co_plus3_dry.pdb          
bond TMC.1.CO TMC.2.N
bond TMC.1.CO TMC.2.N1
bond TMC.1.CO TMC.3.N
bond TMC.1.CO TMC.3.N1
bond TMC.1.CO TMC.4.N
bond TMC.1.CO TMC.4.N1
loadamberparams       /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/BMIM.frcmod       
BMI   = loadmol2      /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/BMIM.mol2  
loadamberparams       /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/NTF2.frcmod       
NSC   = loadmol2      /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/NTF2.mol2  
SYS = loadpdb /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/Co_plus3_bpy3_IL_solvated.pdb      
bond SYS.1.CO SYS.2.N
bond SYS.5.CO SYS.6.N
bond SYS.1.CO SYS.2.N1
bond SYS.5.CO SYS.6.N1
bond SYS.1.CO SYS.3.N
bond SYS.5.CO SYS.7.N
bond SYS.1.CO SYS.3.N1
bond SYS.5.CO SYS.7.N1
bond SYS.1.CO SYS.4.N
bond SYS.5.CO SYS.8.N
bond SYS.1.CO SYS.4.N1
bond SYS.5.CO SYS.8.N1
addions SYS Cl- 0 
set SYS box {32.00,32.00,32.00}
check SYS        
savepdb SYS /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/Co_plus3_bpy3_IL_solvated.pdb        
saveamberparm SYS /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/Co_plus3_bpy3_IL_solvated.prmtop /home/fren5/AutoSolvae-update/AutoSolvate/test_tmc/Co_plus3_bpy3_IL_solvated.inpcrd     
quit              
