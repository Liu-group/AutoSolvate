Tutorial
=============================

Following code walkthrough illustrates the usage of Autosolvate.

There will be two full example systems: napthalene in water and xx in chlorform

Prerequisites
-------------------------------------------
In order to follow along with this tutorial, you need:

autosolvate environment, including amber and babel

autosolvate.py

solute coordinates (napthalene and xx)

>>> conda activate autosolvate

>>> vim napthalene.xyz

>>>  18
>>>  Cartesian coordinates for napthalene, optimized b973c
>>>  C   2.43657685723130     -0.69213325436462      0.00002078124010
>>>  C   2.43658018238548      0.69212401386934      0.00002721987684
>>>  C   1.22901070352146      1.39397634349173      0.00001994687146
>>>  C   0.00000343342043      0.71253090912213      0.00000714065790
>>>  C   -1.22900043327981      1.39398218579512      0.00000100629271
>>>  C   -2.43657340957243      0.69213571136675     -0.00001177104717
>>>  C   -2.43657674766550     -0.69212157876293     -0.00001810130367
>>>  C   -1.22900727009014     -1.39397382444367     -0.00001197564250
>>>  C   0.00000003650493     -0.71252846610373      0.00000054622907
>>>  C   1.22900387312110     -1.39397976876823      0.00000729993441
>>>  H   1.23663547013727     -2.47558382748379      0.00000171706349
>>>  H   -1.23664377992167     -2.47557790006419     -0.00001669525228
>>>  H   -3.37040599960840     -1.23506616836952     -0.00002743124154
>>>  H   -3.37040001518104      1.23508486720593     -0.00001671418790
>>>  H   -1.23663195806702      2.47558626549296      0.00000651593895
>>>  H   1.23664734635303      2.47558038342078      0.00002358751257
>>>  H   3.37040941277112      1.23506867161039      0.00003754491782
>>>  H   3.37040339793988     -1.23508256301448      0.00002608213973


Step 1: Solvate system
-------------------------------------------
Bash commands::

>>> python autosolvate.py -m napthalene.xyz -s water -c 0 -u 1 -g "amber" -o nap_in_water

>>> python autosolvate.py -m xx.xyz -s chloroform -c 1 -u 2 -g "gaussian" -o xx_in_chcl3

Flag definitions:

-m main, solute xyz file

-s solvent, name of solvent (water, methanol, chloroform, nma)

-c charge, formal charge of solute (integer)

-u multiplicity, spin multiplicity of solute (integer)

-g charge method, name of charge fitting method (amber, gaussian)

-o output, desired name of solvated system

-b cube size, size of solvent cube in angstroms (number)

-r srun, option to use slurm to submit as a job

** A note on charge fitting methods **

The semi-emperical charge fitting available through Amber performs well for open-shell systems. However, it is not sufficient for closed-shell systems, which will require the use of quantum chemistry charge fitting methods. The methods currently available are bcc fitting in Amber and resp in Gaussian.

Step 2: Equilibrate and generate QM/MM trajectory
-----------------------------------------------------

Bash commands::

>>> python generatetrajs.py


Step 3: Microsolvated cluster extraction
------------------------------

Bash commands::

>>> python extract.py
