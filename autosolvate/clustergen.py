import getopt, sys, os
import subprocess
import mdtraj as md
import numpy as np


def clustergen(filename='water_solvated', trajname='water_solvated.netcdf', startframe=0, interval=100, size=2, srun_use=False):
    r"""
    Extract microsolvated cluster around center solute

    Parameters
    ----------
    filename : str, Optional, default: 'water_solvated'
        Filename prefix for .prmtop and output files
    trajname :  str, Optional, default: 'water_solvated.netcdf'
        Name of trajectory
    startframe : int, Optional, default: 0
        First frame to extract the microsolvated clusters from trajectory
    interval : int, Optional, default: 100
        Interval at which to extract microsolvated clusters from trajectory
    size : float, Optional, default: 2
        size of solvent shell around center solute in Angstrom
    srun_use : bool, Optional, default: False
        Run all commands with a srun prefix.

    Returns
    -------
    None
        Result stored as 'filename'-cluster.xyz file
    """
    
    print('Loading trajectory')
    a=md.load(trajname, top=filename+'.prmtop')
    traj2=a.xyz
    if frame>-1 and traj2.shape[0]<frame:
      print("trajectory too short")
    print('selecting center solute')
    molecules=a.topology.find_molecules()
    center_list=[]
    for atom in molecules[0]:
      center_list.append(atom.index)
    print('calculating distance to all solvent molecules')
    center_xyz=traj2[frame, np.array(center_list),:]
    center_xyz=np.average(center_xyz,axis=0)
    unit_cell_list=a.unitcell_lengths[0]
    tshape=traj2.shape
    unit_cell_list3 = np.repeat(unit_cell_list[np.newaxis, :], tshape[1], axis=0)
    shift=-center_xyz+unit_cell_list3/2.
    traj3=np.remainder(traj2[frame,:,:]+shift,unit_cell_list3)
    dist_molecules=np.empty(len(molecules))
    select_molecules=[]
    size_molecules=[]
    use_molecules=0
    for i in range(len(molecules)):
      select_true=False
      dist_atom=10000
      for atom in molecules[i]:
        for center_idx in center_list:
            dist4=np.linalg.norm(traj3[atom.index,:]-traj3[center_idx,:])
            if dist4<dist_atom:
                dist_atom=dist4
      dist_molecules[i]=dist_atom
    print('select solvent molecules')
    dist_use=np.sort(dist_molecules)
    ncutout=np.argwhere(dist_use>size)[0][0]
    if ncutout>1:
      ncutout=ncutout+1
    print("selected", ncutout, 'solvent molecules')
    select_mol=np.sort(np.argsort(dist_molecules)[:ncutout])
    print('saving xyz')
    select_list=[]
    for i in select_mol:
      for atom in molecules[i]:
        select_list.append(atom.index)
    select_list.sort()
    print("select atoms", select_list)
    b=a.slice(frame)
    select_xyz=traj3[select_list,:]
    c=b.atom_slice(select_list)
    c.xyz[0,:,:]=select_xyz
    c.save_xyz(filename+'-cutoutn.xyz')

def startclustergen(argumentList):
    print(argumentList)
    options = "f:t:a:i:s:r"
    long_options = ["filename", "trajname", "startframe", "interval", "size", "srunuse"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    srun_use=False
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-f", "-filename"):
            print ("Filename:", currentValue)
            filename=str(currentValue)
        elif currentArgument in ("-t", "-trajname"):
            print ("Trajectory name:", currentValue)
            trajname=str(currentValue)
        elif currentArgument in ("-a", "-startframe"):
            print ("startframe to extract:", currentValue)
            frame=int(currentValue)
        elif currentArgument in ("-i", "-interval"):
            print ("interval to extract:", currentValue)
            interval=int(currentValue)
        elif currentArgument in ("-s", "-size"):
            print ("Cutout size in Angstrom:", currentValue)
            size=float(currentValue)
        elif currentArgument in ("-r", "-srunuse"):
            print("usign srun")
            srun_use=True

    clustergen(filename=filename, trajname=trajname, startframe=startframe, interval=interval, size=size, srun_use=srun_use)

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startclustergen(argumentList)
