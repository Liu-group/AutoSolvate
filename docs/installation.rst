Installation
=============================
AutoSolvate can be installed with pip. 

Dependencies
-----------------

Ensure you have pip and conda installed. Install dependencies openbabel, packmol and Ambertools:

   >>> conda install -c conda-forge openbabel
   >>> conda install -c conda-forge packmol
   >>> conda install -c conda-forge ambertools=21 compilers 

Alternative ways to install Ambertools is described `here <https://ambermd.org/AmberTools.php>`_.

Pip install
----------------

Install autosolvate as following from the commandline:

>>> pip install autosolvate

From source
---------------
Alternative to pip install is installing from source. Download the source code from github.::

   >>> git clone git@github.com:Liu-group/AutoSolvate.git

Go inside the AutoSolvate directory and install it:: 


   python setup.py install


Check
----------------

Check your python installation. These commands in python shouldn't give any errors.::

   from openbabel import pybel
   from openbabel import openbabel as ob
   import autosolvate


Check your ambertools and packmol installation as well::

   >>> which packmol
   >>> which tleap

