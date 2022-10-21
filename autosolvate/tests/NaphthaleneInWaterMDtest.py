## the production files are from the tutorial https://github.com/Liu-group/AutoSolvate/blob/cb04b4a75366ab9103d5aa92accb587f89bfcfe3/docs/tutorial.rst
## autosolvate boxgen -m naphthalene_neutral.xyz -s water -c 0 -u 1 -g "bcc" -o nap_neutral
## autosolvate mdrun -f water_solvated -q 0 -u 1 -t 300 -p 1 -m 10000 -n 10000 -o 100 -s 100 -l 250 -r
import mdtraj as md

def test_Autosolvate_mmmin():

    InfileData = open('mmmin.in','r').readlines()
    OutfileData = open('mmmin.out','r').readlines()
    for linecount, line in enumerate(OutfileData):
        if 'FINAL RESULTS' in line:
            count = linecount + 5
    NSTEP_out =  int(OutfileData[count].split()[0])

    for line in InfileData:
        if 'maxcyc' in line:
            Max_NSTEP_in = int(line.split(',')[0].split('=')[1])
   # print (NSTEP_out, Max_NSTEP_in)
    
    assert NSTEP_out <= Max_NSTEP_in

def test_Autosolvate_mmheat():
    InfileData = open('mmheat.in','r').readlines()
    InfofileData = open('mmheat.info','r').readlines()
    Traj = md.load('water_solvated-heat.netcdf',top='water_solvated.prmtop')
   # print('length of traj',len(Traj))


    for line in InfileData:
        if 'temp0'in line:
           # print(line.split(',')[0].split('=')[1])
            Temp_in = float(line.split(',')[0].split('=')[1])
        elif 'nstlim' in line:
            Nstlim = int(line.split(',')[0].split('=')[1])
        elif 'dt' in line:
            dt = float(line.split(',')[0].split('=')[1])
        elif 'ntpr' in line:
            ntpr = int(line.split(',')[0].split('=')[1])
    
    Time_infile = dt*Nstlim
    Time_netcdf = dt*ntpr*len(Traj)

    for line in InfofileData:
        if 'TIME' in line:
            if 'PRESS' in line:
                Time_outinfo = float(line.split()[5])
                Temp_outinfo = float(line.split()[8])
        
    
    assert (Temp_outinfo - Temp_in) <= 5
    assert Time_infile == Time_outinfo
    assert Time_netcdf == Time_infile

    return (Time_infile)

def test_Autosolvate_mmnpt():
    Time_start = test_Autosolvate_mmheat()
    InfileData = open('mmnpt.in','r').readlines()
    InfofileData = open('mmnpt.info','r').readlines()

    Traj = md.load('water_solvated-mmnpt.netcdf',top='water_solvated.prmtop')

    for line in InfileData:
        if 'temp0'in line:
          #  print(line.split(',')[0].split('=')[1])
            Temp_in = float(line.split(',')[0].split('=')[1])
        elif 'nstlim' in line:
            Nstlim = int(line.split(',')[0].split('=')[1])
        elif 'dt' in line:
            dt = float(line.split(',')[0].split('=')[1])
        elif 'ntpr' in line:
            ntpr = int(line.split(',')[0].split('=')[1])
    
    Time_netcdf = dt*ntpr*len(Traj)

    for line in InfofileData:
        if 'TIME' in line:
            if 'PRESS' in line:
                Time_outinfo = float(line.split()[5]) - Time_start
                Temp_outinfo = float(line.split()[8])
   
    Time_infile = dt*Nstlim
  #  print(Time_infile)
    assert (Temp_outinfo - Temp_in) <= 5
    assert Time_infile == Time_outinfo
    assert Time_netcdf == Time_infile



test_Autosolvate_mmmin()
test_Autosolvate_mmheat()
test_Autosolvate_mmnpt()