Installation
=============================
AutoSolvate can be installed with conda. 

Dependencies
-----------------

* If you install AutoSolvate :ref:`from conda <condainstall>`, you don't need to worry about the dependencies of AutoSolvate because they are automatically installed.

* If you choose to install AutoSolvated :ref:`from source code <sourceinstall>`, you can still install all dependencies automatically, which we will explain :ref:`later in this document <sourceinstall>`.

However, if you are curious about the dependencies of AutoSolvate, please take a look at ``devtools/conda-envs/test_env.yaml`` in the AutoSolvate source code directory. This `YAML <https://yaml.org/>`_ file summarizes all dependencies of AutoSolvate. These packages will allow you to use most functionalities of AutoSolvate.

However, if you'd like to use AutoSolvate to directly drive QM/MM calculations, you need to install TeraChem. `TeraChem <http://www.petachem.com/>`_ is a comercial quantum chemistry software, so you need to purchase and install separately. 

AutoSolvate use the AmberTools to run classical MD simulations without GPU acceleration. If you'd like to use the GPU accelerated version of Amber, please refer to `Amber website <https://ambermd.org/AmberTools.php>`_.

.. _condainstall:

Conda install
----------------

Install autosolvate as following from the commandline::

   >>> conda install -c conda-forge autosolvate

.. _sourceinstall:

From source
---------------
Alternative to conda install is installing from source. Download the source code from github.::

   >>> git clone git@github.com:Liu-group/AutoSolvate.git

First set up a conda environment for AutoSolvate with the needed dependencies::

   >>> conda env create -f devtools/conda-envs/test_env.yaml

This will automatically create a conda environment called ``autosolvate``. Now activate this environment::

   >>> conda activate autosolvate

Go inside the AutoSolvate directory and install it:: 

   >>> python setup.py install


Check
----------------

Check your python installation. These commands in python shouldn't give any errors.::

   from openbabel import pybel
   from openbabel import openbabel as ob
   import autosolvate


Check your ambertools and packmol installation as well::

   >>> which packmol
   >>> which tleap

