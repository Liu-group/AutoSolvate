Tutorial
=============================

Following code walkthrough illustrates the usage of Autosolvate.

There will be two full example systems: napthalene in water and napthalene radical in chloroform

Prerequisites
-------------------------------------------
Once you have AutoSolvate and all dependencies installed you will need the solute xyz file and then you are ready to go! Make sure to give each molecule its own directory to avoid the possibility of overwriting the amber files when running two at a time. The napthalene neutral and radical coordinates are provided below so that you can follow along on your own computer:

::

  18
  napthalene neutral 
       C     2.4397703245   -0.7099883961    0.0000206200
       C     2.4397218526    0.7099981201    0.0000271508
       C     1.2475921776    1.4061556571    0.0000203110
       C    -0.0000128759    0.7189947033    0.0000073141
       C    -1.2476290200    1.4061688746    0.0000008829
       C    -2.4397553974    0.7100487925   -0.0000117263
       C    -2.4397460082   -0.7099448889   -0.0000182422
       C    -1.2476288777   -1.4062156405   -0.0000121401
       C     0.0000138676   -0.7190995078    0.0000006641
       C     1.2476602178   -1.4062240260    0.0000074983
       H     1.2448250471   -2.4927306634    0.0000020169
       H    -1.2447711187   -2.4927196649   -0.0000168971
       H    -3.3840069825   -1.2452230520   -0.0000277743
       H    -3.3839437525    1.2454155894   -0.0000167697
       H    -1.2448430780    2.4926825384    0.0000062499
       H     1.2447883528    2.4926610011    0.0000242506
       H     3.3839630326    1.2452901872    0.0000373621
       H     3.3840333383   -1.2452476243    0.0000259290

::

  18
  napthalene radical
       C     2.4584929186   -0.6980401434    0.0000208854
       C     2.4584830542    0.6980208281    0.0000273558
       C     1.2392834454    1.4064616303    0.0000201346
       C    -0.0000127820    0.7187236077    0.0000072068
       C    -1.2393189424    1.4064428097    0.0000009527
       C    -2.4585398474    0.6980627613   -0.0000119130
       C    -2.4584830245   -0.6980052848   -0.0000182990
       C    -1.2392726206   -1.4064393494   -0.0000121035
       C     0.0000166810   -0.7186826023    0.0000005670
       C     1.2392855074   -1.4064696461    0.0000073561
       H     1.2470800358   -2.4919577916    0.0000018617
       H    -1.2470920207   -2.4919275393   -0.0000168975
       H    -3.3951566422   -1.2429180456   -0.0000277271
       H    -3.3952681112    1.2428765068   -0.0000168560
       H    -1.2469606339    2.4919363439    0.0000063915
       H     1.2471333000    2.4919523490    0.0000239494
       H     3.3951743890    1.2429028846    0.0000376679
       H     3.3951863936   -1.2429173191    0.0000261673

Now that you have the structures, make a neutral and radical directory. We will start with the neutral molecule. 

Step 1: Solvate system
-------------------------------------------

The first step is putting the solute in the solvent box, which uses the boxgen command. The documentation shows all of the options for this command, but the only one that is required is specifying the solute xyz file. It will be listed as -m for main. To run boxgen with all of the defaults, use the following command:
::

  autosolvate boxgen -m napthalene_neutral.xyz 

Autosolvate will use the default values of water as the solvent, solute charge of 0, solute multiplicity of 1, charge fitting method of resp, box size of 54, and output file name of water_solvated. 

If AutoSolvate is running successfully, the following messages will be printed to your screen:
::

  AutoSolvate is starting in command line mode!
  Running the module to generate solvent box and force field parameters.
  ['-m', 'nap_neutral.xyz']
  Main/solutexyz nap_neutral.xyz
  WARNING: Amber home directory is not specified in input options
  WARNING: Checking AMBERHOME environment virable...
  ['echo', '$AMBERHOME']
  WARNING: AMBERHOME detected:  $AMBERHOME
  
  Converting xyz to pdb
  Generate frcmod file for the solute
  cleaning up solute.xyz.pdb
  Then write out mol2
  
  Welcome to antechamber 21.0: molecular input file processor.
  
  acdoctor mode is on: check and diagnose problems in the input file.
  The atom type is set to gaff; the options available to the -at flag are
      gaff, gaff2, amber, bcc, and sybyl.
  -- Check Format for pdb File --
     Status: pass
  -- Check Unusual Elements --
     Status: pass
  -- Check Open Valences --
     Status: pass
  -- Check Geometry --
       for those bonded   
       for those not bonded   
     Status: pass
  -- Check Weird Bonds --
     Status: pass
  -- Check Number of Units --
     Status: pass
  acdoctor mode has completed checking the input file.
  
  Info: Total number of electrons: 68; net charge: 0
  
  Running: /jet/home/agale/miniconda3/envs/autosolvate/bin/sqm -O -i sqm.in -o sqm.out
  
  Finally generate frcmod with parmchk2
  Now create the solute library file
  Generate Amber parameters for the solvated system
  Now add pre-equlibrated solvent box to the solute
  The script has finished successfully

Additionally, you should now have the following files in your directory:
::

  ANTECHAMBER_AC.AC           ATOMTYPE.INF              nap_neutral.xyz   sqm.in   
  ANTECHAMBER_AC.AC0          leap_add_solventbox.cmd   solute.frcmod     sqm.out  
  ANTECHAMBER_AM1BCC.AC       leap_add_solventbox.log   solute.lib        sqm.pdb  
  ANTECHAMBER_AM1BCC_PRE.AC   leap.cmd                  solute.mol2       water_solvated.inpcrd
  ANTECHAMBER_BOND_TYPE.AC    leap.log                  solute.pdb        water_solvated.pdb
  ANTECHAMBER_BOND_TYPE.AC0   leap_savelib.log          solute.xyz.pdb    water_solvated.prmtop

The three files that we care about for moving forward to the next step are water_solvated.pdb, water_solvated.inpcrd, and water_solvated.prmtop.

The .inpcrd file contains the input coordinates, and the .prmtop file contains the Amber paramter topology. The .pdb file has the coordinates for the solute in the solvent box, so you want to check that both the solvent and the solute are there:
::

        CRYST1   66.461   66.696   66.822  90.00  90.00  90.00 P 1           1
        ATOM      1  C   SLU     1       2.302  -0.634   0.016  1.00  0.00
        ATOM      2  C1  SLU     1       2.302   0.786   0.016  1.00  0.00
        ATOM      3  C2  SLU     1       1.110   1.482   0.016  1.00  0.00
        ATOM      4  C3  SLU     1      -0.138   0.795   0.016  1.00  0.00
        ATOM      5  C4  SLU     1      -1.386   1.482   0.016  1.00  0.00
        ATOM      6  C5  SLU     1      -2.578   0.786   0.016  1.00  0.00
        ATOM      7  C6  SLU     1      -2.578  -0.634   0.016  1.00  0.00
        ATOM      8  C7  SLU     1      -1.386  -1.330   0.016  1.00  0.00
        ATOM      9  C8  SLU     1      -0.138  -0.643   0.016  1.00  0.00
        ATOM     10  C9  SLU     1       1.110  -1.330   0.016  1.00  0.00
        ATOM     11  H   SLU     1       1.107  -2.417   0.016  1.00  0.00
        ATOM     12  H1  SLU     1      -1.383  -2.417   0.016  1.00  0.00
        ATOM     13  H2  SLU     1      -3.522  -1.169   0.016  1.00  0.00
        ATOM     14  H3  SLU     1      -3.522   1.321   0.016  1.00  0.00
        ATOM     15  H4  SLU     1      -1.383   2.569   0.016  1.00  0.00
        ATOM     16  H5  SLU     1       1.107   2.569   0.016  1.00  0.00
        ATOM     17  H6  SLU     1       3.246   1.321   0.016  1.00  0.00
        ATOM     18  H7  SLU     1       3.246  -1.169   0.016  1.00  0.00
        TER
        ATOM     19  O   WAT     2      30.753  27.440  26.571  1.00  0.00
        ATOM     20  H1  WAT     2      30.672  26.525  26.300  1.00  0.00
        ATOM     21  H2  WAT     2      30.339  27.937  25.865  1.00  0.00
        TER
        ATOM     22  O   WAT     3      28.885  29.218  28.452  1.00  0.00
        ATOM     23  H1  WAT     3      28.109  28.738  28.742  1.00  0.00
        ATOM     24  H2  WAT     3      29.536  28.538  28.277  1.00  0.00

The fourth column has 18 'SLU' entries, or solvent, and under that there are 6 'WAT' entries, which we can see makes up two water molecules. 

With these three files, we are ready to proceed to the next step!




In order to change some of these or make sure everything is defined explicitly, we can use more of the flag options.

>>> python autosolvate.py -m napthalene_neutral.xyz -s water -c 0 -u 1 -g "bcc" -o nap_neutral_water

Now, the charge fitting has been changed to semi-emperical bcc fitting in Amber, which is appropriate for a closed-shell system, and the output file name is more specific. The solvent, charge, and multiplicity have all been defined explicitly as well.



** A note on charge fitting methods **

The semi-emperical charge fitting available through Amber performs well for closed-shell systems. However, it is not sufficient for open-shell systems, which will require the use of quantum chemistry charge fitting methods. The methods currently available are bcc fitting in Amber and resp in Gaussian.

Step 2: Equilibrate and generate QM/MM trajectory
-----------------------------------------------------

The second step is running QM/MM, which includes equilibration and production time. For this tutorial, we will run a very fast demonstration just to see how the mdrun command works.

To do a short example run of QM/MM use the following command:
::

  autosolvate mdrun -f water_solvated -q 0 -u 1 -t 300 -p 1 -m 10000 -n 10000 -o 100 -s 100 -l 10 -r "True"
  
The mdrun command has several more options than the previous one, but the only required options are filename, charge, and multiplicity (the first three in the command above). It is also important to use the srun option if you are using HPC or having a queueing system, because otherwise it will run on the login node and will exceed the acceptable time limit.

If AutoSolvate is running successfully, the following messages will be printed to your screen:
::

  AutoSolvate is starting in command line mode!
  Running the module to automatically run MD simulations of solvated structure.
  ['-f', 'water_solvated', '-q', '0', '-u', '1', '-t', '300', '-p', '1', '-m', '10000', '-n', '10000', '-o', '100', '-s', '100', '-l', '10', '-r', 'True']
  Filename: water_solvated
  Charge: 0
  Spinmultiplicity: 1
  Temperature in K: 300
  Pressure in bar: 1
  Steps MM heat: 10000
  Steps MM NPT: 10000
  Steps QMMM heat: 100
  Steps QMMM NPT: 100
  Steps QMMM min: 10
  using srun
  MM Energy minimization
  srun: job 5791719 queued and waiting for resources
  srun: job 5791719 has been allocated resources
  MM Heating
  srun: job 5791725 queued and waiting for resources
  srun: job 5791725 has been allocated resources



The main output here is the QM/MM trajectory nap_neutral_water-qmmmnvt.netcdf.




Longer MM and QM/MM steps are necessary to reach equilibration, and the default settings are more appropriate than what is used here for a production run. The default mdrun will have the following settings:

MM heat
    temperature=300 K
    stepsmmheat=10000 steps
MM NPT
    pressure=1 bar
    stepsmmnpt=300000 steps
QM/MM
    charge=0
    spinmult=1
QM/MM min
    stepsqmmmmin=250 steps
QM/MM heat
    stepsqmmmheat=1000 steps
QM/MM NVT
    stepsqmmmnvt=10000 steps
    
When you are ready to do a production run and want to use all of these defaults, you can use the dry run option to generate the input files without running them to make sure that everything looks right: 
::

  autosolvate mdrun -f water_solvated -q 0 -u 1 -d
  
If AutoSolvate is running correctly, it will print the following messages:
::

  AutoSolvate is starting in command line mode!
  Running the module to automatically run MD simulations of solvated structure.
  ['-f', 'water_solvated', '-q', '0', '-u', '1', '-d']
  Filename: water_solvated
  Charge: 0
  Spinmultiplicity: 1
  Dry run mode: only generate the commands to run MD programs and save them into a file without executing the commands
  MM Energy minimization
  MM Heating
  MM NPT equilibration
  QMMM Energy minimization
  QMMM Heating
  QMMM NVT Run
  
The following files will be added to your directory:
::

  mmheat.in  qmmmheat.in  runMM.sh
  mmmin.in   qmmmmin.in   runQMMMM.sh
  mmnpt.in   qmmmnvt.in   tc_job.tpl



Step 3: Microsolvated cluster extraction
----------------------------------------------------------

Bash commands to extract 4 Angstrom solvent shell for each 10th frame or every 5fs:

>>> autosolvate clustergen -f nap_neutral_water -t nap_neutral_water-qmmmnvt.netcdf -a 0 -i 10 -s 4

Main output are the microsolvated clusters ``nap_neutral_water-cutoutn-*.xyz``.


Second System: Napthalene Radical
----------------------------------------------------------

>>> autosolvate.py boxgen -m napthalene_radical.xyz -s chloroform -c 1 -u 2 -g "resp" -o nap_radical_chcl3



