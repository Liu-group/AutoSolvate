Tutorial
=============================

Following code walkthrough illustrates the usage of Autosolvate.

Step 1: Solvate system
-------------------------------------------
Bash commands::

>>> python autosolvate.py -m solute.xyz -s water


Step 2: Equilibrate and generate QM/MM trajectory
-----------------------------------------------------

Bash commands::

>>> python generatetrajs.py


Step 3: Microsolvated cluster extraction
------------------------------

Bash commands::

>>> python extract.py
