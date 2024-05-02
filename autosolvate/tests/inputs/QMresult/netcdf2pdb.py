import mdtraj as md

# 加载轨迹和拓扑文件
traj = md.load('Fe_plus2_solvated-mmnpt.netcdf', top='Fe_plus2_solvated.prmtop')

# 保存轨迹的第一帧为PDB文件
traj[-1].save_pdb('Fe_plus2_solvated-mmnpt-last.pdb')
