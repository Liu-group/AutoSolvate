import getopt, sys, os
import subprocess




def runMM(filename='water_solvated', minin='min.in', heatin='heat.in', mmin='mm.in', srun_use=False):
    r"""
    Equilibrate with MM
   
    Parameters
    ----------
    filename : str, Optional, default: 'water_solvated'
        Filename prefix for .prmtop and inpcrd files
    minin : str, Optional, default: 'min.in'
        Amber input file for energy minimization
    minin : str, Optional, default: 'min.in'
        Amber input file for energy minimization
    heatin : str, Optional, default: 'heat.in'
        Amber input file for energy minimization
    mmin : str, Optional, default: 'mm.in'
        Amber input file for MM equilibration
    srun_use : bool, Optional, default: False
        Run all commands with a srun prefix.
    
    Returns
    -------
    None
        Results stored in .netcdf files and log files
    """
    print('MM Energy minimzation')
    cmd='sander -O -i '+minin+' -o min.out -p '+filename+'.prmtop -c '+filename+'.inpcrd -r min.ncrst'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('MM Heating')
    cmd='sander -O -i '+heatin+' -o heat.out -p '+filename+'.prmtop -c min.ncrst -r heat.ncrst -x '+filename+'-heat.netcdf -inf heat.mdinfo'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('MM equilibration')
    cmd='sander -O -i '+mmin+' -o md.out -p '+filename+'.prmtop -c heat.ncrst -r md.ncrst -x '+filename+'.netcdf -inf md.info'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)


def runQMMM(filename='water_solvated', qmmmminin='qmmmmin.in', qmmmrunin='qmmmrun.in', spin_mult=1, srun_use=False):
    r"""
    Equilibrate and then generate trajectory with QM/MM

    Parameters
    ----------
    None
    
    Returns
    -------
    None
        Results stored in .netcdf files and log files
    """
    cmd=""
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)



if __name__ == '__main__':
    argumentList = sys.argv[1:]
    print(argumentList)
    options = "m:f:h:d:u:r"
    arguments, values = getopt.getopt(argumentList, options, long_options)
    srun_use=False
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-f", "-filename"):
            print ("Filename:", currentValue)
            filename=str(currentValue)
        elif currentArgument in ("-m", "-minin"):
            print ("Minin:", currentValue)
            minin=str(currentValue)
        elif currentArgument in ("-h", "-heatin"):
            print ("Heatin:", currentValue)
            heatin=str(currentValue)
        elif currentArgument in ("-d", "-mmin"):
            print ("MMin:", currentValue)
            mmin=str(currentValue)
        elif currentArgument in ("-q", "-qmmmminin"):
            print ("QM/MM Min in:", currentValue)
            mdin=str(currentValue)
        elif currentArgument in ("-n", "-qmmmrunin"):
            print ("QM/MM Run in:", currentValue)
            mdin=str(currentValue)
        elif currentArgument in ("-u", "-spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            spinmult=int(currentValue)
        elif currentArgument in ("-r", "-srunuse"):
            print("usign srun")
            srun_use=True

     
    runMM(filename=filename, minin=minin, heatin=heatin, mmin=mmin, srun_use=srun_use)
    runQMMM(filename=filename, qmmmminin=minin, qmmmrunin=qmmmrunin, spinmult=spinmult, srun_use=srun_use)
