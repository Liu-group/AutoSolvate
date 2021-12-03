import getopt, sys, os
import subprocess


def writeminin():
    r"""
    Write Amber MM minimization input file 
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
        Result stored as min.in
    """
    f = open("min.in","w")
    f.write("Minimize\n")
    f.write("&cntrl\n")
    f.write("imin=1,\n")
    f.write("ntx=1,\n")
    f.write("maxcyc=2000,\n")
    f.write("ncyc=1000,\n")
    f.write("ntpr=100,\n")
    f.write("ntwx=0,\n")
    f.write("cut=8.0,\n")
    f.write("/\n")
    f.close()


def writeheatin(temperature=298.15, stepsheat=10000):
        r"""
        Write Amber MM heat input file 

        Parameters
        ----------
        temperature : float, Optional, default: 298.15
            temperature to heat to
        stepsheat : int, Optional, default: 10000
            MM steps for heating.

        Returns
        -------
        None
            Result stored as heat.in
        """
        f = open("heat.in","w")
        f.write("Heat\n")
        f.write("&cntrl\n")
        f.write("imin=0,\n")
        f.write("ntx=1,\n")
        f.write("nstlim="+str(stepsheat)+",\n")
        f.write("dt=0.002,\n")
        f.write("ntf=2,\n")
        f.write("ntc=2,\n")
        f.write("tempi=0.0,\n")
        f.write("temp0="+str(temperature)+",\n")
        f.write("ntpr=100,\n")
        f.write("ntwx=100,\n")
        f.write("cut=8.0,\n")
        f.write("ntb=1,\n")
        f.write("ntp=0,\n")
        f.write("ntt=3,\n")
        f.write("gamma_ln=2.0,\n")
        f.write("nmropt=1,\n")
        f.write("ig=-1,\n")
        f.write("/\n")
        f.write("&wt type='TEMP0', istep1=0, istep2="+str(stepsheat)+", value1=0.0, value2="+str(temperature)+" /\n")
        f.write("&wt type='END' /\n")
        f.close()

def runMM(filename='water_solvated', srun_use=False):
    r"""
    Equilibrate with MM
   
    Parameters
    ----------
    filename : str, Optional, default: 'water_solvated'
        Filename prefix for .prmtop and inpcrd files
    srun_use : bool, Optional, default: False
        Run all commands with a srun prefix.
    
    Returns
    -------
    None
        Results stored as .netcdf files and log files
    """
    print('MM Energy minimization')
    cmd='sander -O -i min.in -o min.out -p '+filename+'.prmtop -c '+filename+'.inpcrd -r min.ncrst'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('MM Heating')
    cmd='sander -O -i heat.in -o heat.out -p '+filename+'.prmtop -c min.ncrst -r heat.ncrst -x '+filename+'-heat.netcdf -inf heat.mdinfo'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('MM equilibration')
    cmd='sander -O -i mm.in -o md.out -p '+filename+'.prmtop -c heat.ncrst -r md.ncrst -x '+filename+'.netcdf -inf md.info'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)

def writeqmmmheatin(temperature=298.15, stepsqmmmheat=1000):
        r"""
        Write QMMM heating input file

        Parameters
        ----------
        temperature : float, Optional, default: 298.15
            temperature to heat to
        stepsqmmmheat : int, Optional, default: 10000
            QMMM steps for heating.

        Returns
        -------
        None
            Result stored qmmmheat.in
        """
        f = open("qmmmheat.in","w")
        f.write("gpr QMMM heat\n")
        f.write(" &cntrl\n")
        f.write("  imin   = 0,\n")
        f.write("  irest  = 1, ! 0- new simulation 1- restart\n")
        f.write("  ntx    = 5, ! 1-read in coordinates, but not velocity, 5-both\n")
        f.write("  cut    = 8.0,\n")
        f.write("  ig     = -1, !random seed\n")
        f.write("  ntc    = 2, ntf    = 2, !Shake is used for solvent\n")
        f.write("  ntb    = 1,\n")
        f.write("  tempi  = "+str(temperature)+",\n")
        f.write("  temp0  = "+str(temperature)+",\n")
        f.write("  ntt    = 3, ! 1=rescale 3=Langevin dynamics 7-bussi\n")
        f.write("  vrand  = 1,\n")
        f.write("  tautp  = 0.01,\n")
        f.write("  gamma_ln = 5.0, !Langevin dynamics collision frequency\n")
        f.write("  nstlim = "+str(stepsqmmmheat)+",  !num steps\n")
        f.write("  dt = 0.0005,  !in ps\n")
        f.write("  ntpr   = 1, !print detials to log every step\n")
        f.write("  ntwx   = 1, !write coordinates to mdcrd every step\n")
        f.write("  ntwr   = 1, !write restart file every step\n")
        f.write("/\n")
        f.write(" &qmmm\n")
        f.write("  qmmask    = ':1',\n")
        f.write("  qmcharge  = 0,\n")
        f.write("  qmshake  = 0, !Shake QM H atoms if shake is turned on (NTC>1) (default)\n")
        f.write("  qm_theory = 'EXTERN',\n")
        f.write("  qm_ewald  = 0,\n")
        f.write("  qmgb      = 0,\n")
        f.write("  verbosity = 2,\n")
        f.write("  writepdb =1,\n")
        f.write(" /\n")
        f.write(" &tc\n")
        f.write("  executable = 'terachem',\n")
        f.write("  use_template = 1,\n")
        f.write(" /\n")
        f.close()


def runQMMM(filename='water_solvated', srun_use=False, spinmult=1):
    r"""
    Run QMMM minimisation, heating, trajectory run

    Parameters
    ----------
    filename : str, Optional, default: 'water_solvated'
        Filename prefix for .prmtop and inpcrd files
    srun_use : bool, Optional, default: False
        Run all commands with a srun prefix.

    Returns
    -------
    None
        Results stored in .netcdf files and log files
    """
    print('QMMM Energy minimization')
    cmd='sander -O -i qmmmmin.in -o qmmmmin.out -p '+filename+'.prmtop -c md.ncrst -r qmmmmin.ncrst  -inf qmmmmin.info -x '+filename+'-min.netcdf'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('For higher Spin multiplicity adjust input file')
    if spinmult!=1:
      cmd="sed -i '3 a guess        ./scr/ca0 ./scr/cb0' tc_job.tpl"
      if self.srun_use:
        cmd='srun -n 1 '+cmd
      subprocess.call(cmd, shell=True)
    cmd='sander -O -i qmmmmin.in -o qmmmmin.out -p '+filename+'.prmtop -c qmmmmin.ncrst -r qmmmmin.ncrst  -inf qmmmmin.info -x '+filename+'-min.netcdf'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('QMMM Heating')
    cmd='sander -O -i qmmmheat.in -o qmmmheat.out -p '+filename+'.prmtop -c qmmmmin.ncrst -r qmmmheat.ncrst  -inf qmmmheat.info -x '+filename+'-heat.netcdf'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    print('QMMM Run')
    cmd='sander -O -i qmmmrun.in -o qmmmrun.out -p '+filename+'.prmtop -c qmmmheat.ncrst -r qmmmrun.ncrst  -inf qmmmrun.info -x '+filename+'-run.netcdf'
    if self.srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)


if __name__ == '__main__':
    argumentList = sys.argv[1:]
    print(argumentList)
    options = "f:t:p:h:s:u:r"
    long_options = ["filename", "temp", "pressure", "stepsheat", "stepsmm", "stepsqmmmheat", "stepsqmmm", "spinmultiplicity", "srunuse"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    srun_use=False
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-f", "-filename"):
            print ("Filename:", currentValue)
            filename=str(currentValue)
        elif currentArgument in ("-t", "-temp"):
            print ("Temperature:", currentValue)
            temperature=float(currentValue)
        elif currentArgument in ("-p", "-pressure"):
            print ("Pressure:", currentValue)
            pressure=float(currentValue)
        elif currentArgument in ("-h","-stepsheat"):
            print ("Steps heat:", currentValue)
            stepsheat=int(currentValue)
        elif currentArgument in ("-s", "-stepsmm"):
            print ("Steps MM:", currentValue)
            stepsmm=int(currentValue)
        elif currentArgument in ("-h","-stepsqmmmheat"):
            print ("Steps QMMM heat:", currentValue)
            stepsqmmmheat=int(currentValue)
        elif currentArgument in ("-s", "-stepsqmmm"):
            print ("Steps QMMM:", currentValue)
            stepsmm=int(currentValue)
        elif currentArgument in ("-u", "-spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            spinmult=int(currentValue)
        elif currentArgument in ("-r", "-srunuse"):
            print("usign srun")
            srun_use=True

    writeminin()
    writeheatin(temperature=temperature, stepsheat=stepsheat)
    writemmin(temperature=temperature, pressure=pressure, stepsmm=stepsmm)
    runMM(filename=filename, srun_use=srun_use)
    writeqmmmminin(spinmult=spinmult)
    writeqmmmheatin(temperature=temperature, stepsqmmmheat=stepsqmmmheat)
    writeqmmmrunin(temperature=temperature, stepsqmmm=stepsqmmm)
    runQMMM(filename=filename, srun_use=srun_use)
