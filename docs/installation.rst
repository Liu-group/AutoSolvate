Installation
=============================
AutoSolvate can be installed with :ref:`from source code <sourceinstall>` or :ref:`conda <condainstall>`

For either approach, you want to make sure to have `git <https://git-scm.com/>`_ and `conda <https://docs.conda.io/en/latest/>`_ installed on your computer.

Dependencies
-----------------

.. note::

   **Windows users**: **AmberTools** cannot be installed from conda, and therefore cannot be automatically installed with the following approaches. **AmberTools** is the only dependency for Windows Users and needs separate installation.

   **Mac/Linux users**: No dependency as long as you follow the instructions below.


* If you install AutoSolvate :ref:`from source code <sourceinstall>`, you can install all dependencies automatically, which we will explain :ref:`later in this document <sourceinstall>`.

* If you choose to install AutoSolvate :ref:`from conda <condainstall>`, you don't need to worry about the dependencies of AutoSolvate because they are automatically installed.

However, if you are curious about the dependencies of AutoSolvate, please take a look at ``devtools/conda-envs/test_env.yaml`` in the AutoSolvate source code directory. This `YAML <https://yaml.org/>`_ file summarizes all dependencies of AutoSolvate. These packages will allow you to use most functionalities of AutoSolvate.

However, to use all functionalities of AutoSolvate, one needs to install a few other packages not included in AutoSolvate installation:

#. If you want to solvate open-shell molecules (spin-multiplicity > 1), RESP charge fitting is needed, which uses `Gaussian <https://gaussian.com/>`_. `Gaussian <https://gaussian.com/>`_ is a commercial quantum chemistry package, so you need to purchase and install separately. 

#. If you'd like to use AutoSolvate to directly drive QM/MM calculations, you need to install TeraChem. `TeraChem <http://www.petachem.com/>`_ is a commercial quantum chemistry package, so you need to purchase and install separately. 

#. AutoSolvate uses the AmberTools to run classical MD simulations without GPU acceleration. If you'd like to use the GPU accelerated version of Amber, please refer to `Amber website <https://ambermd.org/AmberTools.php>`_.


.. _sourceinstall:

From source
---------------
We recommend installing AutoSolvate from source code because of its good compatability cross operating systems, except for some Windows, which needs some extra care (more info :ref:`below <windowsinstallwarning>`).

Download the source code from github.::

   >>> git clone git@github.com:Liu-group/AutoSolvate.git

First set up a conda environment for AutoSolvate with the needed dependencies::

   >>> conda env create -f devtools/conda-envs/test_env.yaml

This will automatically create a conda environment called ``autosolvate``. Now activate this environment::

   >>> conda activate autosolvate

Go inside the AutoSolvate directory and install it:: 

   >>> python setup.py install

   
.. _windowsinstallwarning:

.. warning::
  
    **Window users**: Since conda installation of AmberTools is not directly
    available on Windnows, you will need to install AmberTools sperately
    from source following the instructions on 
    `Amber Documentation <https://ambermd.org/GetAmber.php#ambertools>`_.

    Alternatively, we may install Windows 10/11 "subsystem for Linux" (WSL/WSL2)
    as instructed on this `Amber webpage <https://ambermd.org/InstWindows.php>`_.
    Once you have WSL/WSL2, you should be able to directly install AutoSolvate
    (including AmberTools dependencies) with the procedure described above for Linux systems.
   
.. _condainstall:   

Conda install
----------------

.. note::

   There are known issues about conda installation of AutoSolvate on **Mac with M1 chips** and **Windows**, and one needs to take extra steps after doing the conda installation. Therefore, you are recommended to install :ref:`from source<sourceinstall>`

Alternative to installing from source is conda install. This works for **Linux** or **old Mac without M1 chips**. Install autosolvate as following from the commandline::

   >>> conda install -c liugroupemory autosolvate

To check out more about the AutoSolvate conda package, please visit `this page on Anaconda.org <https://anaconda.org/LiuGroupEmory/autosolvate>`_.

Following are workarounds for **Mac with M1 chips** and **Windows**.

**Mac with M1 chips**:

To make the conda installation of AutoSolvate work on Mac with the M1 chip, one needs to take care of two things:

#. Download **mini-forge**. Then do ``conda install -c liugroupemory autosolvate``. See more info `<https://stackoverflow.com/questions/65534042/anaconda-and-upgrading-to-new-m1-mac>`_ or `<https://towardsdatascience.com/using-conda-on-an-m1-mac-b2df5608a141>`_.
#. You need to do conda installation of **AmberTools** again by typing ``conda install -c conda-forge ambertools``. 

**Windows**: 

You can still do conda installation of AutoSolvate on Windows, but AmberTools does not get automatically installed after you run::

>>> conda install -c liugroupemory autosolvate

So you need to install AmberTools separately. You can follow the intrucitons :ref:`here<windowsinstallwarning>`.



Check
----------------

Check your python installation. These commands in python shouldn't give any errors::

   from openbabel import pybel
   from openbabel import openbabel as ob
   import autosolvate


Check your ambertools and packmol installation as well::

   >>> which packmol
   >>> which tleap

