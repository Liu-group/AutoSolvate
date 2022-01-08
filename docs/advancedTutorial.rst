Advanced Tutorial
=============================
Here we introduce some advanced usages of AutoSolvate. To learn the basic usages, please refer to the basic :doc:`tutorial` page.

Advanced Example 1: Custom Solvent
------------------------------------
Apart from the 5 common solvents contained in AutoSolvate, the user can use custom solvents found in databases or literature to build the solvated structure, as long as the `.frcmod` and `.off` files are available.

Here we show a simple use case. We are still going to work on the neutral napthalene molecule used in the basic :doc:`tutorial`. However, this time we will put it in a custom solvent, dimethylsulfoxide (DMSO), which is not contained in AutoSolvate.


Step 1: Find custom solvent force field files
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Theoretically, we can generate GAFF force field for any solvent, and use that for our simulation.

However, it is ideal to use the solvent force fields that have been verified in publications and can reliably reproduce experimental results.

For example, if you want to simulate some solute in DMSO, you may want to look for existing Amber force field of DMSO. One online resource is the AMBER parameter database hosted by the Bryce Group at the University of Manchester: `<http://amber.manchester.ac.uk/>`_.

On the website, you may find some solvent boxes available, including DMSO. You will find two downloadable files for the DMSO solvent box:

#. **OFF**: The DMSO solvent box library file 
#. **FRCMOD**: The DMSO force field modification file

You can download them and save as: `dmso.frcmod` and `dmso.off`.

Next, it is very important to find out the name of the solvent box file that will be recognized by AmberTools, and pass this name to AutoSolvate.

If you open the `dmso.off` file, you will see the first few lines as below:

.. code-block:: text
   :linenos:

   !!index array str
    "d"
   !entry.d.unit.atoms table  str name  str type  int typex  int resx  int flags  int seq  int elmnt  dbl chg
    "S" "S" 0 1 131073 1 16 0.307524
    "CT1" "CT" 0 1 131073 2 6 -0.262450


Notice that for a solvent box `OFF` file, the 2nd line is the **name** of the solvent box that can be recognized by AmberTool/tleap.
In this case, the name is `d`. That means, if one loads this `OFF` file and uses tleap to add the solvent box, the corresponding command should be::

>>> solventbox [solute_unit_name] d [box_size] [closeness]

So `d` is the solvent name that we should pass to AutoSolvate.

.. note::

   If you don't like the original solvent box name given in the OFF file, feel free to change it to something else. For example, you can change the 2nd line
   of dmso.off to "DMSO". Later you want to pass the new name, "DMSO" to AutoSolvate


Step 2: Run AutoSolvate with the custom solvent
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

To generate the solvent box structure and MD prmtop files with custom solvent, the basic procedure is the same as the simple exmaple about adding water (see :doc:`tutorial`).

The only difference is to provide 3 extra options:
#. solvent name (not the real name, but the name given in `OFF` file) with option `-s`
#. solvent `OFF` file path with option `-l`
#. solvent `FRCMOD` file path with option `-p`

Assuming that you have `napthalene_neutral.xyz`, `dmso.off`, `dmso.frcmod` files all in the current
working directory, and the environment with AutoSolvate installed has been actived.
To add the DMSO solvent box to the neutral napthalene molecule, you can simply run the following command::

>>> autosolvate boxgen -m napthalene_neutral.xyz -s d  -l dmsobox.off -p dmso.frcmod

This command should generate the solvated files: `d_solvated.inpcrd`, `d_solvated.prmtop`, and `d_solvated.pdb`
