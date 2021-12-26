import getopt, sys, os
import subprocess


def writeMMminInput():
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
    f.write("maxcyc=2000,\n")
    f.write("ncyc=1000,\n")
    f.write("ntpr=100,\n")
    f.write("ntwx=0,\n")
    f.write("cut=8.0,\n")
    f.write("/\n")
    f.close()


def writeMMheatInput(temperature=300, stepsmmheat=10000):
        r"""
        Write Amber MM heat input file 

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
        f.write("nmropt=1,\n")
        f.write("ig=-1,\n")
        f.write("/\n")
        f.write("&wt type='TEMP0', istep1=0, istep2="+str(stepsmmheat)+", value1=0.0, value2="+str(temperature)+" /\n")
        f.write("&wt type='END' /\n")
        f.close()

def writeMMNPTInput(temperature=300, pressure=1, stepsmmnpt=300000):
        r"""
        Write Amber MM NPT input file

        Parameters
        ----------
        temperature : float, Optional, default: 300
            temperature 
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
        f.write("imin=0,\n")
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
        f.write("pres0 = "+str(pressure)+"\n,")
        f.write("ntt=3,\n")
        f.write("gamma_ln=2.0,\n")
        f.write("ig=-1,\n")
        f.write("/\n")
        f.close()



def runMM(filename='water_solvated', stepsmmheat=10000, stepsmmnpt=300000, stepsmmnvt=0, srun_use=False, pmemduse=False):
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
    cmd='sander -O -i mmmin.in -o mmmin.out -p '+filename+'.prmtop -c '+filename+'.inpcrd -r mmmin.ncrst -inf mmmin.info'
    if srun_use:
      cmd='srun -n 1 '+cmd
    subprocess.call(cmd, shell=True)
    if stepsmmheat>0:
      print('MM Heating')
      cmd='sander -O -i mmheat.in -o mmheat.out -p '+filename+'.prmtop -c mmmin.ncrst -r mmheat.ncrst -x '+filename+'-heat.netcdf -inf mmheat.info'
      if srun_use:
        cmd='srun -n 1 '+cmd
      subprocess.call(cmd, shell=True)
    if stepsmmnpt>0:
      print('MM NPT equilibration')
      cmd=' -O -i mmnpt.in -o mmnpt.out -p '+filename+'.prmtop -c mmheat.ncrst -r mm.ncrst -x '+filename+'-mmnpt.netcdf -inf mmnpt.info'
      if pmemduse:
        cmd= 'pmemd.cuda' +cmd
      else:
        cmd= 'sander'+ cmd
      if srun_use:
        cmd='srun -n 1 '+cmd
      subprocess.call(cmd, shell=True)


def writeQMMMTemplate(spinmult=1,charge=0):
        r"""
        Write Terachem file-based interface to Amber

        Parameters
        ----------
        spinmult : float, Optional, default: 1
            spin multiplicity of solvated system
        charge : int, Optional, default: 0
            charge of solvated system

        Returns
        -------
        None
            Result stored as mmnpt.in
        """
        f = open("tc_job.tpl","w")
        f.write("basis        lacvps_ecp\n")
        if spinmult==1:
          f.write("method       b3lyp\n")
        else: 
          f.write("method       ub3lyp\n")
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
        temperature : float, Optional, default: 300
            temperature to heat to
        stepsqmmmmin : int, Optional, default: 250
            QMMM steps for minimise.

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

def writeQMMMInput(temperature=300, charge=0, stepsqmmm=250, infilename='qmmmheat.in'):
        r"""
        Write QMMM heating input file

        Parameters
        ----------
        temperature : float, Optional, default: 300
            temperature to heat to
        stepsqmmm : int, Optional, default: 250
            QMMM steps 

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
        f.write("  ntt    = 3, ! 1=rescale 3=Langevin dynamics 7-bussi\n")
        f.write("  vrand  = 1,\n")
        f.write("  tautp  = 0.01,\n")
        f.write("  gamma_ln = 5.0, !Langevin dynamics collision frequency\n")
        f.write("  nstlim = "+str(stepsqmmm)+",  !num steps\n")
        f.write("  dt     = 0.0005,  !in ps\n")
        f.write("  ntpr   = 1, !print detials to log every step\n")
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


def runQMMM(filename='water_solvated', spinmult=1, srun_use=False, stepsqmmmmin=250, stepsqmmmheat=1000, stepsqmmmnvt=10000):
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
    if stepsqmmmmin>0:
      print('QMMM Energy minimization')
      cmd='sander -O -i qmmmmin.in -o qmmmmin.out -p '+filename+'.prmtop -c mm.ncrst -r qmmmmin.ncrst  -inf qmmmmin.info -x '+filename+'-qmmmmin.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      subprocess.call(cmd, shell=True)
      if spinmult>1:
        print('Adjusting terachem input file for higher Spin multiplicity')
        cmd="sed -i '3 a guess        ./scr/ca0 ./scr/cb0' tc_job.tpl"
        if srun_use:
          cmd='srun -n 1 '+cmd
        subprocess.call(cmd, shell=True)
        cmd='sander -O -i qmmmmin.in -o qmmmmin.out -p '+filename+'.prmtop -c qmmmmin.ncrst -r qmmmmin.ncrst  -inf qmmmmin.info -x '+filename+'-qmmmmin.netcdf'
        if srun_use:
          cmd='srun -n 1 '+cmd
        subprocess.call(cmd, shell=True)
    if stepsqmmmheat>0:
      print('QMMM Heating')
      cmd='sander -O -i qmmmheat.in -o qmmmheat.out -p '+filename+'.prmtop -c qmmmmin.ncrst -r qmmmheat.ncrst  -inf qmmmheat.info -x '+filename+'-qmmmheat.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      subprocess.call(cmd, shell=True)
    if stepsqmmmnvt>0:
      print('QMMM NVT Run')
      cmd='sander -O -i qmmmnvt.in -o qmmmnvt.out -p '+filename+'.prmtop -c qmmmheat.ncrst -r qmmmnvt.ncrst  -inf qmmmnvt.info -x '+filename+'-qmmmnvt.netcdf'
      if srun_use:
        cmd='srun -n 1 '+cmd
      subprocess.call(cmd, shell=True)


def startmd(argumentList):
    print(argumentList)
    options = "f:t:p:m:n:l:o:s:q:u:r:x"
    long_options = ["filename", "temp", "pressure", "stepsmmheat", "stepsmmnpt", "stepsqmmmmin", "stepsqmmmheat", "stepsqmmmnvt", "charge", "spinmultiplicity", "srunuse", "pmemduse"]
    arguments, values = getopt.getopt(argumentList, options, long_options)
    srun_use=False
    temperature=300
    pressure=1
    stepsmmheat=10000
    stepsmmnpt=300000
    stepsqmmmmin=250
    stepsqmmmheat=1000
    stepsqmmmnvt=10000
    pmemduse=False
    for currentArgument, currentValue in arguments:
        if currentArgument in ("-f", "-filename"):
            print ("Filename:", currentValue)
            filename=str(currentValue)
        elif currentArgument in ("-t", "-temp"):
            print ("Temperature in K:", currentValue)
            temperature=float(currentValue)
        elif currentArgument in ("-p", "-pressure"):
            print ("Pressure in bar:", currentValue)
            pressure=float(currentValue)
        elif currentArgument in ("-m","-stepsmmheat"):
            print ("Steps MM heat:", currentValue)
            stepsmmheat=int(currentValue)
        elif currentArgument in ("-n", "-stepsmmnpt"):
            print ("Steps MM NPT:", currentValue)
            stepsmmnpt=int(currentValue)
        elif currentArgument in ("-l","-stepsqmmmmin"):
            print ("Steps QMMM min:", currentValue)
            stepsqmmmmin=int(currentValue)
        elif currentArgument in ("-o","-stepsqmmmheat"):
            print ("Steps QMMM heat:", currentValue)
            stepsqmmmheat=int(currentValue)
        elif currentArgument in ("-s", "-stepsqmmmnvt"):
            print ("Steps QMMM NPT:", currentValue)
            stepsqmmmnvt=int(currentValue)
        elif currentArgument in ("-q", "-charge"):
            print("Charge:", currentValue)
            charge=int(currentValue)
        elif currentArgument in ("-u", "-spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            spinmult=int(currentValue)
        elif currentArgument in ("-r", "-srunuse"):
            print("usign srun")
            srun_use=True
        elif currentArgument in ("-x", "-pmemduse"):
            print("usign pmemd.cuda instead of sander")
            srun_use=True


    writeMMminInput()
    writeMMheatInput(temperature=temperature, stepsmmheat=stepsmmheat)
    writeMMNPTInput(temperature=temperature, pressure=pressure, stepsmmnpt=stepsmmnpt)
    
    writeQMMMMinInput(stepsqmmmmin=stepsqmmmmin)
    writeQMMMTemplate(spinmult=spinmult, charge=charge)
    writeQMMMInput(temperature=temperature, stepsqmmm=stepsqmmmheat, charge=charge, infilename='qmmmheat.in')
    writeQMMMInput(temperature=temperature, stepsqmmm=stepsqmmmnvt, charge=charge, infilename='qmmmnvt.in' )
    
    runMM(filename=filename, stepsmmheat=stepsmmheat, stepsmmnpt=stepsmmnpt, srun_use=srun_use, pmemduse=pmemduse)
    
    runQMMM(filename=filename, spinmult=spinmult, srun_use=srun_use, stepsqmmmmin=stepsqmmmmin, stepsqmmmheat=stepsqmmmheat, stepsqmmmnvt=stepsqmmmnvt)

if __name__ == '__main__':
    argumentList = sys.argv[1:]
    startmd(argumentList)
