import getopt, sys, os
import subprocess



def writeminin():
    r"""
    Write amber minimization input file 

    Parameters
    ----------
    None

    Returns
    -------
    None
        Result stored min.in
    """
        f = open("min.in","w")
        f.write("Minimize\n")
        f.write("&cntrl\n")
        f.write("imin=1,\n")
        f.write("ntx=1,\n")
        f.write("maxcyc=2000,\n")
        f.write("ncyc=1000,\n")
        f.write("set SLU tail SLU.1."+self.tail + "\n")
        f.write("ntpr=100,\n")
        f.write("ntwx=0,\n")
        f.write("cut=8.0,\n")
        f.write("/\n")
        f.close()

def writeheatin(temperature=298.15, stepsheat=10000):
    r"""
    Write amber heat input file 

    Parameters
    ----------
    temperature : float, Optional, default: 298.15
        temperature to heat to
    stepsheat : int, Optional, default: 10000
        MM steps for heating.

    Returns
    -------
    None
        Result stored heat.in
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
        Results stored in .netcdf files and log files
    """
    print('MM Energy minimzation')
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
    options = "f:t:p:h:s:u:r"
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
        elif currentArgument in ("-n", "-qmmmrunin"):
            print ("QM/MM Run in:", currentValue)
            mdin=str(currentValue)
        elif currentArgument in ("-u", "-spinmultiplicity"):
            print ("Spinmultiplicity:", currentValue)
            spinmult=int(currentValue)
        elif currentArgument in ("-r", "-srunuse"):
            print("usign srun")
            srun_use=True

    writeminin()
    writeheatin(temperature=temperature, stepsheat=stepsheat)
    writemmin(temperature=temperature, pressure=pressure, stepsmm=stepsmm)
    runMM(filename=filename, minin=minin, heatin=heatin, mmin=mmin, srun_use=srun_use)
    runQMMM(filename=filename, qmmmminin=minin, qmmmrunin=qmmmrunin, spinmult=spinmult, srun_use=srun_use)
