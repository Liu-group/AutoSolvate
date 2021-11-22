.. autosolvate documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AutoSolvate's documentation!
=====================================================================
*Automated workflow to solvate molecules and run QM/MM trajectories*

Description
---------------------

This package enables automated initial structure generation for explicitly solvate systems. This includes input file preparation. Additionally automated QM/MM trajectory generation is supported for the explicitly solvated systems. These features empower the user to rapidly generate large computational data sets.


Dependencies
-------------------

Currently Autosolvate depends on following packages:

#. `Open Babel <http://openbabel.org/>`_
#. `AmberTools21 <https://ambermd.org/AmberTools.php>`_
#. `Packmol <http://leandro.iqm.unicamp.br/m3g/packmol/home.shtml>`_

Optionally `Gaussian <http://gaussian.com/>`_ can be used to estimate partial charges.

For the generation of QM/MM following packages are supported:

#. `TeraChem <http://www.petachem.com/products.html>`_

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   api



