import getopt, sys, os
import subprocess


def writeMMminInput(stepsmmmin=2000):
    r"""
    Write Amber MM minimization input file 
    
    Parameters
    ----------
    None
    
    Returns
    -------
    None
        Result stored as mmmin.in
    """
    f = open("mmmin.in","w")
    f.write("Minimize\n")
    f.write("&cntrl\n")
    f.write("imin=1,\n")
    f.write("ntx=1,\n")
    f.write("maxcyc="+str(stepsmmmin)+",\n")
    f.write("ncyc=1000,\n")
    f.write("ntpr=100,\n")
    f.write("ntwx=0,\n")
    f.write("cut=8.0,\n")
    f.write("/\n")
    f.close()


def writeMMheatInput(temperature=300, stepsmmheat=10000):
        r"""
        Write Amber MM heating input file 

        Parameters
        ----------
        temperature : float, Optional, default: 300
            temperature to heat to
        stepsmmheat : int, Optional, default: 10000
            MM steps for heating.

        Returns
        -------
        None
            Result stored as mmheat.in
        """
        f = open("mmheat.in","w")
        f.write("Heat\n")
        f.write("&cntrl\n")
        f.write("imin=0,\n")
        f.write("ntx=1,\n")
        f.write("nstlim="+str(stepsmmheat)+",\n")
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
        f.write("ig=-1,\n")
        f.write("/\n")
        f.write("&wt type='TEMP0', istep1=0, istep2="+str(stepsmmheat)+", value1=0.0, value2="+str(temperature)+" /\n")
        f.write("&wt type='END' /\n")
        f.close()


def writeMMNVEInput(stepsmmnve=10000):
        r"""
        Write Amber MM NVE input file

        Parameters
        ----------
        stepsmmnve : int, Optional, default: 10000
            Number of MM steps for NVE.

        Returns
        -------
        None
            Result stored as mmnve.in
        """
        f = open("mmnve.in","w")
        f.write("NVE\n")
        f.write("&cntrl\n")
        f.write("imin=0, irest=1,\n")
        f.write("ntx=5,\n")
        f.write("nstlim="+str(stepsmmnve)+",\n")
        f.write("dt=0.002,\n")
        f.write("ntf=2,\n")
        f.write("ntc=2,\n")
        f.write("ntpr=100,\n")
        f.write("ntwx=100,\n")
        f.write("cut=8.0,\n")
        f.write("ntb=1,\n")
        f.write("ntp=0,\n")
        f.write("ntt=0,\n")
        f.write("ig=-1,\n")
        f.write("/\n")
        f.close()


def writeMMNPTInput(temperature=300, pressure=1, stepsmmnpt=300000):
        r"""
        Write Amber MM NPT input file

        Parameters
        ----------
        temperature : float, Optional, default: 300
            temperature in Kelvin
        pressure : float, Optional, default: 1
            pressure in bar
        stepsmmnpt : int, Optional, default: 10000
            MM NPT steps

        Returns
        -------
        None
            Result stored as mmnpt.in
        """
        f = open("mmnpt.in","w")
        f.write("MM NPT\n")
        f.write("&cntrl\n")
        f.write("imin=0, irest=1,\n")
        f.write("ntx=5,\n")
        f.write("nstlim="+str(stepsmmnpt)+",\n")
        f.write("dt=0.002,\n")
        f.write("ntf=2,\n")
        f.write("ntc=2,\n")
        f.write("temp0="+str(temperature)+",\n")
        f.write("ntpr=1000,\n")
        f.write("ntwx=1000,\n")
        f.write("cut=8.0,\n")
        f.write("ntb=2,\n")
        f.write("ntp=1,\n")
        f.write("taup=1,\n") 
        f.write("pres0 = "+str(pressure)+",\n")
        f.write("ntt=3,\n")
        f.write("gamma_ln=2.0,\n")
        f.write("ig=-1,\n")
        f.write("/\n")
        f.close()



def runMM(filename='water_solvated', stepsmmheat=10000, stepsmmnve=0, stepsmmnpt=300000, srun_use=False, pmemduse=False, dryrun=False):
    r"""
    Equilibrate with MM
   
    Parameters
    ----------
    filename : str, Optional, default: 'water_solvated'
        Filename prefix for .prmtop and .inpcrd files
    stepsmmheat : int, Optional, default: 10000
        MM heating steps
    stepsmmnve : int, Optional, default: 0
        MM NVE steps
    stepsmmnpt : int, Optional, default: 300000
        MM NPT steps
    srun_use : bool, Optional, default: False
        Run all commands with a srun prefix.
    pmemduse :  bool, Optional, default: False
        Use pmemd.CUDA instead of sander
    dryrun : bool, Optional, default: False
        Dry run mode: only generate the commands to run MD programs and save them into a file without executing the commands
    
    Returns
    -------
    None
        Results stored as .netcdf files and log files
    """
    frun = None
    if dryrun:
       frun = open('runMM.sh','w')
        
    print('MM Energy minimization')
    cmd=' -O -i mmmin.in -o mmmin.out -p '+filename+'.prmtop -c '+filename+'.inpcrd -r mm.ncrst -inf mmmin.info'
    if pmemduse:
        cmd= 'pmemd.cuda' +cmd
    else:
        cmd= 'sander'+ cmd
    if srun_use:
      cmd='srun -n 1 '+cmd

    if dryrun:
        frun.write(cmd+'\n')
    else:
        subprocess.call(cmd, shell=True)
    if stepsmmheat>0:
      print('MM Heating')
      cmd=' -O -i mmheat.in -o mmheat.out -p '+filename+'.prmtop -c mm.ncrst -r mm.ncrst -x '+filename+'-heat.netcdf -inf mmheat.info'
      if pmemduse:
        cmd= 'pmemd.cuda' +cmd
      else:
        cmd= 'sander'+ cmd
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)
    if stepsmmnve>0:
      print('MM NVE equilibration')
      cmd=' -O -i mmnve.in -o mmnve.out -p '+filename+'.prmtop -c mm.ncrst -r mm.ncrst -x '+filename+'-mmnve.netcdf -inf mmnve.info'
      if pmemduse:
        cmd= 'pmemd.cuda' +cmd
      else:
        cmd= 'sander'+ cmd
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)
    if stepsmmnpt>0:
      print('MM NPT equilibration')
      cmd=' -O -i mmnpt.in -o mmnpt.out -p '+filename+'.prmtop -c mm.ncrst -r mm.ncrst -x '+filename+'-mmnpt.netcdf -inf mmnpt.info'
      if pmemduse:
        cmd= 'pmemd.cuda' +cmd
      else:
        cmd= 'sander'+ cmd
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)

    if dryrun:
       frun.close()

def writeQMMMTemplate(spinmult=1,charge=0,functional="b3lyp"):
        r"""
        Write Terachem file-based interface to Amber

        Parameters
        ----------
        spinmult : float, default: 1
            spin multiplicity of solvated system
        charge : int, Optional, default: 0
            charge of solvated system

        Returns
        -------
        None
            Result stored as tc_job.tpl
        """
        f = open("tc_job.tpl","w")
        f.write("basis        lacvps_ecp\n")
        if spinmult==1:
          f.write("method       "+functional+"\n")
        else: 
          f.write("method       u"+functional+"\n")
        f.write("dispersion   yes\n")
        f.write("scf          diis+a\n")
        f.write("threall      1e-13\n")
        f.write("convthre     2e-5\n")
        f.write("xtol         1-4\n")
        f.write("dftd         no\n")
        f.write("maxit        200\n")
        f.write("dftgrid      1\n")
        f.write("charge       "+str(charge)+"\n")
        f.write("spinmult     "+str(spinmult)+"\n")
        f.write("scrdir       ./scr\n")
        f.write("keep_scr     yes\n")
        f.write("end\n")
        f.close()

def writeQMMMMinInput(stepsqmmmmin=250):
        r"""
        Write QMMM min input file

        Parameters
        ----------
        stepsqmmmmin : int, Optional, default: 250
            QMMM steps for minimization

        Returns
        -------
        None
            Result stored qmmmmin.in
        """
        f = open("qmmmmin.in","w")
        f.write("gpr QMMM heat\n")
        f.write(" &cntrl\n")
        f.write("  imin   = 0,\n")
        f.write("  irest  = 0, ! 0- new simulation 1- restart\n")
        f.write("  ntx    = 1, ! 1-read in coordinates, but not velocity, 5-both\n")
        f.write("  cut    = 8.0,\n")
        f.write("  ig     = -1, !random seed\n")
        f.write("  ntc    = 2, ntf    = 2, !Shake is used for solvent\n")
        f.write("  ntb    = 1,\n")
        f.write("  tempi  = 0,\n")
        f.write("  temp0  = 0,\n")
        f.write("  ntt    = 3, ! 1=rescale 3=Langevin dynamics 7-bussi\n")
        f.write("  vrand  = 1,\n")
        f.write("  tautp  = 0.01,\n")
        f.write("  gamma_ln = 500.0, !Langevin dynamics collision frequency\n")
        f.write("  nstlim = "+str(stepsqmmmmin)+",  !num steps\n")
        f.write("  dt = 0.0001,  !in ps\n")
        f.write("  ntpr   = 1, !print detials to log every step\n")
        f.write("  ntwx   = 1, !write coordinates to mdcrd every step\n")
        f.write("  ntwr   = 1, !write restart file every step\n")
        f.write("  ifqnt  = 1,!turn on QM/MM\n")
        f.write("/\n")
        f.write(" &qmmm\n")
        f.write("  qmmask    = ':1',\n")
        f.write("  qmcut     = 8,\n")
        f.write("  qmcharge  = 0,\n")
        f.write("  qmshake   = 0, !Shake QM H atoms if shake is turned on (NTC>1) (default)\n")
        f.write("  qm_theory = 'EXTERN',\n")
        f.write("  qm_ewald  = 0,\n")
        f.write("  qmgb      = 0,\n")
        f.write("  verbosity = 2,\n")
        f.write("  writepdb  = 1,\n")
        f.write(" /\n")
        f.write(" &tc\n")
        f.write("  executable = 'terachem',\n")
        f.write("  use_template = 1,\n")
        f.write(" /\n")
        f.close()

def writeQMMMInput(temperature=300, charge=0, stepsqmmm=250, infilename='qmmmheat.in', nve=False):
        r"""
        Write QMMM heating or NPT trajectory input file

        Parameters
        ----------
        temperature : float, Optional, default: 300
            temperature to heat to
        charge : int, Optional, default: 0
            Total charge of system
        stepsqmmm : int, Optional, default: 250
            QMMM steps 
        infilename : str, Optional, default: 'qmmmheat.in'
            Filename to save Amber input file
        nve : bool, Optional, default: False
            Set True if running nve

        Returns
        -------
        None
            Result stored in infilename
        """
        f = open(infilename,"w")
        f.write("gpr QMMM "+infilename+"\n")
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
        if nve==True:
          f.write("  ntt    = 0, ! 0=const. E, 1=rescale 3=Langevin dynamics 7-bussi\n")
        else:
          f.write("  ntt    = 3, ! 0=const. E, 1=rescale 3=Langevin dynamics 7-bussi\n")
        f.write("  vrand  = 1,\n")
        f.write("  tautp  = 0.01,\n")
        if nve!=True:
          f.write("  gamma_ln = 5.0, !Langevin dynamics collision frequency\n")
        f.write("  nstlim = "+str(stepsqmmm)+",  !num steps\n")
        f.write("  dt     = 0.0005,  !in ps\n")
        f.write("  ntpr   = 1, !print details to log every step\n")
        f.write("  ntwx   = 1, !write coordinates to mdcrd every step\n")
        f.write("  ntwr   = 1, !write restart file every step\n")
        f.write("  ifqnt  = 1,!turn on QM/MM\n")
        f.write("/\n")
        f.write(" &qmmm\n")
        f.write("  qmmask    = ':1',\n")
        f.write("  qmcharge  = "+str(charge)+",\n")
        f.write("  qmcut     = 8,\n")
        f.write("  qmshake   = 0, !Shake QM H atoms if shake is turned on (NTC>1) (default)\n")
        f.write("  qm_theory = 'EXTERN',\n")
        f.write("  qm_ewald  = 0,\n")
        f.write("  qmgb      = 0,\n")
        f.write("  verbosity = 2,\n")
        f.write("  writepdb  = 1,\n")
        f.write(" /\n")
        f.write(" &tc\n")
        f.write("  executable = 'terachem',\n")
        f.write("  use_template = 1,\n")
        f.write(" /\n")
        f.close()


def runQMMM(filename='water_solvated', spinmult=1, srun_use=False, stepsqmmmmin=250, stepsqmmmheat=1000, stepsqmmmnve=0, stepsqmmmnvt=10000, dryrun=False):
    r"""
    Run QMMM minimization, heating and NVT trajectory run

    Parameters
    ----------
    filename : str, Optional, default: 'water_solvated'
        Filename prefix for .prmtop input file and .netcdf output file
    spinmult : int, Required
        Spin multiplicity of system
    srun_use : bool, Optional, default: False
        Run all commands with a srun prefix.
    stepsqmmmmin : int, Optional, default: 250
        Number of QMMM minimization steps
    stepsqmmmheat : int, Optional, default: 1000
        Number of QMMM heating steps
    stepsqmmmnve : int, Optional, default: 0
        Number of QMMM NVE steps
    stepsqmmmnvt : int, Optional, default: 10000
        Number of QMMM NVT trajectory steps
    dryrun : bool, Optional, default: False
        Dry run mode: only generate the commands to run MD programs and save them into a file without executing the commands

    Returns
    -------
    None
        Results stored in .netcdf files and log files
    """
    frun = None
    if dryrun:
       frun = open('runQMMMM.sh','w')

    if stepsqmmmmin>0:
      print('QMMM Energy minimization')
      cmd='sander -O -i qmmmmin.in -o qmmmmin.out -p '+filename+'.prmtop -c mm.ncrst -r qmmm.ncrst  -inf qmmmmin.info -x '+filename+'-qmmmmin.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)
      if spinmult>1:
        print('Adjusting terachem input file for higher Spin multiplicity')
        cmd="sed -i '3 a guess        ./scr/ca0 ./scr/cb0' tc_job.tpl"
        if srun_use:
          cmd='srun -n 1 '+cmd
        if dryrun:
          frun.write(cmd+'\n')
        else:
          subprocess.call(cmd, shell=True)
        cmd='sander -O -i qmmmmin.in -o qmmmmin.out -p '+filename+'.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmmin.info -x '+filename+'-qmmmmin.netcdf'
        if srun_use:
          cmd='srun -n 1 '+cmd
        if dryrun:
          frun.write(cmd+'\n')
        else:
          subprocess.call(cmd, shell=True)
    if stepsqmmmheat>0:
      print('QMMM Heating')
      cmd='sander -O -i qmmmheat.in -o qmmmheat.out -p '+filename+'.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmheat.info -x '+filename+'-qmmmheat.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)
    if stepsqmmmnve>0:
      print('QMMM NVE Run')
      cmd='sander -O -i qmmmnve.in -o qmmmnve.out -p '+filename+'.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmnve.info -x '+filename+'-qmmmnve.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)
    if stepsqmmmnvt>0:
      print('QMMM NVT Run')
      cmd='sander -O -i qmmmnvt.in -o qmmmnvt.out -p '+filename+'.prmtop -c qmmm.ncrst -r qmmm.ncrst  -inf qmmmnvt.info -x '+filename+'-qmmmnvt.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      if dryrun:
        frun.write(cmd+'\n')
      else:
        subprocess.call(cmd, shell=True)

    if dryrun:
        frun.close()


def startmd(argumentList):
    r"""
    Wrap function that parses command line options for autosolvate clustergen,
    generates inputfiles for Amber and TeraChem,
    runs MM and QMMM stages.

    Parameters
    ----------
    argumentList: list
        The list contains the command line options to specify MM and QMMM stage options.

        Command line options:
          -f, --filename  prefix of .prmtop and .inpcrd files
          -t, --temp  temperature in Kelvin to equilibrate
          -p, --pressure  pressure in bar to equilibrate during MM NPT step
          -i, --stepsmmmin Number of MM minimization steps
          -m, --stepsmmheat  Number of MM heating steps, setting to 0 skips the MM heating step
          -b, --stepsmmnve  Number of MM NVE steps, setting to 0 skips the MM NVE step
          -n, --stepsmmnpt  Number of MM NPT steps, setting to 0 skips the MM NPT step
          -l, --stepsqmmmmin  Number of QMMM minimization steps, setting to 0 skips the QMMM minimization step
          -o, --stepsqmmmheat  Number of QMMM heating steps, setting to 0 skips the QMMM heating step
          -v, --stepsqmmmnve  Number of QMMM NVE steps, setting to 0 skips the QMMM NVE step
          -s, --stepsqmmmnvt  Number of QMMM NVT steps, setting to 0 skips the QMMM NVT step
          -q, --charge  Total charge of system
          -u, --spinmultiplicity  Spin multiplicity of whole system
          -k, --functional  DFT functional to use for the QM part in QM/MM
          -r, --srunuse  option to run inside a slurm job
          -x, --pmemduse  Speed up MM with pmemd.CUDA instead of sander
          -d, --dryrun  Dry run mode: only generate the commands to run MD programs and save them into a file without executing the command
          -h, --help  short usage description    

    Returns
    -------
    None
        Generate MD simulation input files and execute MD programs, or same the MD program execution commands in runMM.sh and runQMMM.sh
    
    Currently some simulation parameters like simulation time step, integrator type, nonbonded cutoff, thermostat type, Langevin collision frequency, barostat type, pressure relaxation time and frequency of trajectory writing can not be changed from default values by the user.

    """
    #print(argumentList)
    options = "hf:t:p:i:m:b:n:l:o:v:s:q:u:k:rxd"
    long_options = ["help", "filename", "temp", "pressure", "stepsmmmin", "stepsmmheat", "stepsmmnve", "stepsmmnpt", "stepsqmmmmin", "stepsqmmmheat", "stepsqmmmnve", "stepsqmmmnvt", "charge", "spinmultiplicity","functional", "srunuse", "pmemduse","dryrun"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    srun_use=False
    temperature=300
    pressure=1
    stepsmmheat=10000
    stepsmmnve=0
    stepsmmnpt=300000
    stepsqmmmmin=250
    stepsqmmmheat=1000
    stepsqmmmnve=0
    stepsqmmmnvt=10000
    charge=0
    spinmult=1
    functional="b3lyp"
    srun_use=False
    pmemduse=False
    dryrun=False
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-h", "-help"):
            print('Usage: autosolvate mdrun [OPTIONS]')
            print('  -f, --filename             prefix of .prmtop and .inpcrd files')
            print('  -t, --temp                 temperature in Kelvin to equilibrate')
            print('  -p, --pressure             pressure in bar to equilibrate during MM NPT step')
            print('  -i, --stepsmmmin           Number of MM minimization steps')
            print('  -m, --stepsmmheat          Number of MM heating steps, setting to 0 skips the MM heating step')
            print('  -b, --stepsmmnve           Number of MM NVE steps, setting to 0 skips the MM NVE step')
            print('  -n, --stepsmmnpt           Number of MM NPT steps, setting to 0 skips the MM NPT step')
            print('  -l, --stepsqmmmmin         Number of QMMM minimization steps, setting to 0 skips the QMMM minimization step')
            print('  -o, --stepsqmmmheat        Number of QMMM heating steps, setting to 0 skips the QMMM heating step')
            print('  -v, --stepsqmmmnve         Number of QMMM NVE steps, setting to 0 skips the QMMM NVE step')
            print('  -s, --stepsqmmmnvt         Number of QMMM NVT steps, setting to 0 skips the QMMM NVT step')
            print('  -q, --charge               Total charge of system')
            print('  -u, --spinmultiplicity     Spin multiplicity of whole system')
            print('  -k, --functional           DFT functional to use for the QM part in QM/MM')
            print('  -r, --srunuse              option to run inside a slurm job')
            print('  -x, --pmemduse             Speed up MM with pmemd.CUDA instead of sander')
            print('  -d, --dryrun               Dry run mode')
            print('  -h, --help                 short usage description')
            exit()
        elif currentArgument in ("-f", "-filename"):
            print ("Filename:", currentValue)
            filename=str(currentValue)
        elif currentArgument in ("-t", "-temp"):
            print ("Temperature in K:", currentValue)
            temperature=float(currentValue)
        elif currentArgument in ("-p", "-pressure"):
            print ("Pressure in bar:", currentValue)
            pressure=float(currentValue)
        elif currentArgument in ("-i","-stepsmmmin"):
            print ("Steps MM min:", currentValue)
            stepsmmmin=int(currentValue)
        elif currentArgument in ("-m","-stepsmmheat"):
            print ("Steps MM heat:", currentValue)
            stepsmmheat=int(currentValue)
        elif currentArgument in ("-b","-stepsmmnve"):
            print ("Steps MM NVE:", currentValue)
            stepsmmnve=int(currentValue)
        elif currentArgument in ("-n", "-stepsmmnpt"):
            print ("Steps MM NPT:", currentValue)
            stepsmmnpt=int(currentValue)
        elif currentArgument in ("-l","-stepsqmmmmin"):
            print ("Steps QMMM min:", currentValue)
            stepsqmmmmin=int(currentValue)
        elif currentArgument in ("-o","-stepsqmmmheat"):
            print ("Steps QMMM heat:", currentValue)
            stepsqmmmheat=int(currentValue)
        elif currentArgument in ("-v","-stepsqmmmnve"):
            print ("Steps QMMM NVE:", currentValue)
            stepsqmmmnve=int(currentValue)
        elif currentArgument in ("-s", "-stepsqmmmnvt"):
            print ("Steps QMMM NPT:", currentValue)
            stepsqmmmnvt=int(currentValue)
        elif currentArgument in ("-q", "-charge"):
            print("Charge:", currentValue)
            charge=int(currentValue)
        elif currentArgument in ("-u", "-spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            spinmult=int(currentValue)
        elif currentArgument in ("-k", "-functional"):
            print ("DFT functional:", currentValue)
            functional=currentValue
        elif currentArgument in ("-r", "-srunuse"):
            print("using srun")
            srun_use=True
        elif currentArgument in ("-x", "-pmemduse"):
            print("using pmemd.cuda instead of sander")
            pmemduse=True
        elif currentArgument in ("-d", "-dryrun"):
            print("Dry run mode: only generate the commands to run MD programs and save them into a file without executing the commands")
            dryrun=True


    writeMMminInput(stepsmmmin=stepsmmmin)
    writeMMheatInput(temperature=temperature, stepsmmheat=stepsmmheat)
    writeMMNVEInput(stepsmmnve=stepsmmnve)
    writeMMNPTInput(temperature=temperature, pressure=pressure, stepsmmnpt=stepsmmnpt)
    
    writeQMMMMinInput(stepsqmmmmin=stepsqmmmmin)
    writeQMMMTemplate(spinmult=spinmult, charge=charge, functional=functional)
    writeQMMMInput(temperature=temperature, stepsqmmm=stepsqmmmheat, charge=charge, infilename='qmmmheat.in')
    writeQMMMInput(stepsqmmm=stepsqmmmnve, charge=charge, infilename='qmmmnve.in', nve=True)
    writeQMMMInput(temperature=temperature, stepsqmmm=stepsqmmmnvt, charge=charge, infilename='qmmmnvt.in' )
    
    runMM(filename=filename, stepsmmheat=stepsmmheat, stepsmmnpt=stepsmmnpt, stepsmmnve=stepsmmnve, srun_use=srun_use, pmemduse=pmemduse, dryrun=dryrun)
    
    runQMMM(filename=filename, spinmult=spinmult, srun_use=srun_use, stepsqmmmmin=stepsqmmmmin, stepsqmmmheat=stepsqmmmheat, stepsqmmmnve=stepsqmmmnve, stepsqmmmnvt=stepsqmmmnvt, dryrun=dryrun)

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startmd(argumentList)
