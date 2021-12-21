Tutorial
=============================

Following code walkthrough illustrates the usage of Autosolvate.

There will be two full example systems: napthalene in water and napthalene radical in chlorform

Prerequisites
-------------------------------------------
In order to follow along with this tutorial, you need:


>>> conda activate autosolvate

>>> vim napthalene_neutral.xyz

::
18
-3.8580765971568428e+02 neutral napthalene
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

>>> vim napthalene_radical.xyz

::
18
-3.8552686324110755e+02 napthalene radical
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
::


Step 1: Solvate system
-------------------------------------------
Bash commands::

>>> python autosolvate.py -m napthalene_neutral.xyz 

Using this command, Autosolvate will initial using the default values of water as the solvent, solute charge of 0, solute multiplicity of 0, charge fitting method of resp, box size of 54, and output file name of water_solvated. In order to change some of these or make sure everything is defined explicitly, we can use more of the flag options.

>>> python autosolvate.py -m napthalene_neutral.xyz -s water -c 0 -u 1 -g "bcc" -o nap_neutral_water

Now, the charge fitting has been changed to semi-emperical bcc fitting in Amber, which is appropriate for a closed-shell system, and the output file name is more specific. The solvent, charge, and multiplicity have all been defined explicitly as well.



** A note on charge fitting methods **

The semi-emperical charge fitting available through Amber performs well for closed-shell systems. However, it is not sufficient for open-shell systems, which will require the use of quantum chemistry charge fitting methods. The methods currently available are bcc fitting in Amber and resp in Gaussian.

Step 2: Equilibrate and generate QM/MM trajectory
-----------------------------------------------------

Bash commands::

>>> python generatetrajs.py


Step 3: Microsolvated cluster extraction
----------------------------------------------------------

Bash commands::

>>> python extract.py


Second System: Napthalene Radical
----------------------------------------------------------

>>> python autosolvate.py -m napthalene_radical.xyz -s chloroform -c 1 -u 2 -g "resp" -o nap_radical_chcl3
>>> python generatetrajs.py
>>> python extract.py



