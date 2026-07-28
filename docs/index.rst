.. autosolvate documentation master file, created by
   sphinx-quickstart on Thu Mar 15 13:55:56 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to AutoSolvate's documentation!
=====================================================================
*Automated workflow to solvate molecules and run QM/MM trajectories*

.. image:: ../autosolvate/GUI/images/logo.png
   :align: center


Description
---------------------

.. note::
   
   |documentation_version_message|


This `open-source package <https://github.com/Liu-group/AutoSolvate/>`_ enables automated initial structure generation for explicitly solvated systems. This includes input file preparation. Additionally automated QM/MM trajectory generation and microsolvated cluster extraction is supported for the explicitly solvated systems. These features empower the user to rapidly generate large computational data sets.

AutoSolvate provides six complementary interfaces that support different user communities and workflows. These interfaces fall into three categories:

**Local Interfaces**

* **JSON input interface** – the recommended interface for constructing complex systems, including organometallic compounds, mixed-solvent systems, and high-throughput workflows.
* **Interactive CLI** – guides users through the setup process with an interactive question-and-answer interface.
* **Legacy command-line interface (CLI)** – maintains compatibility with AutoSolvate 1.0 command-line workflows.
* **Legacy GUI** – preserves the complete graphical interface from AutoSolvate 1.0 for backward compatibility and existing workflows.

**AI Interface**

* **MCP interface** – enables large language model (LLM) agents to interact with AutoSolvate through the Model Context Protocol (MCP).

**Cloud Interface**

* **AutoSolvateWeb** – a cloud-based, chatbot-assisted platform that enables users to build solvated systems through natural-language conversations without local software installation.

Quickstart
----------------------
#. To start using AutoSolvate, you may read the :doc:`installation` page.
#. To run simple jobs with a single organic molecule as the solute species and a single solvent species, you may start with the following options:
        ##. To use the command line interface (CLI), you may read the :doc:`tutorial`.
        ##. To use the graphical user interface (GUI), you may read the :doc:`tutorialGUI`.
#. To generate solvated organometallic compounds, you may read the :doc:`tutorial_for_metalcomplexes`.
#. To generate solvated metal ions, you may read :ref:`prepare-solvated-ion`.
#. To build a mixed solvent system, you may read the :doc:`tutorialMulticomponent`.
#. To use the AutoSolvateWeb, you may read the :doc:`tutorialWeb` and :doc:`webJobParameters`.
#. To use AutoSolvate as an MCP for AI workflows, you may read the :doc:`mcp_setup`
#. To import AutoSolvate as a Python API, you may refer to the :doc:`api`.
 

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   installation
   tutorial
   tutorialWeb
   tutorialGUI
   tutorialMulticomponent
   tutorial_for_metalcomplexes
   mcp_setup
   webJobParameters
   api
   citation



