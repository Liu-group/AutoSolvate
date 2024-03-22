AutoSolvateWeb Job Parameters
=============================

Available job parameters for the web interface are listed in the following Table.

.. list-table:: **General Parameter**
   :widths: auto
   :header-rows: 1
   :class: longtable

   * - **Parameter**
     - **Description**
     - **Type**
     - **Default Value**
     - **Valid Range**
   * - Solute
     - The molecule or group of molecules of interest. The file path to solute xyz file is accepted.
     - str
     - ""
     - Not given
   * - Solvent
     - The substance used to mimic the environment surrounding a solute molecule. Available solvents: 'water', 'methanol', 'chloroform', 'nma', 'acetonitrile'.
     - str
     - 'water'
     - NA
   * - Charge Method
     - Method that determines partial atomic charges. Available methods: 'bcc', 'resp'. Use 'resp' (quantum mechanical calculation needed) or 'bcc' to estimate partial charges.
     - str
     - 'resp'
     - 'resp' only for open-shell system
   * - Solute Charge
     - Net charge of solute, the solvent box will be neutralized with Cl- and Na+ ions.
     - int
     - 0
     - NA
   * - Solute Spin Multiplicity
     - Spin multiplicity of solute. Defined as 2S+1, where S is the total spin quantum number of the molecule. A practical guide about calculating S for molecules is available [here](https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Supplemental_Modules_(Physical_and_Theoretical_Chemistry)/Electronic_Structure_of_Atoms_and_Molecules/Evaluating_Spin_Multiplicity).
     - int
     - 1
     - S >= 1
   * - Solute Cube Size
     - Size of the solvent box for molecular dynamics simulations (in Angstrom).
     - int
     - 54
     - L >=20
   * - Dry Run
     - Only generate the input files and batch job scripts to run MD programs without executing.
     - Bool
     - False
     - NA
   * - Temperature
     - Temperature in Kelvin to equilibrate in MM or QM/MM calculation.
     - int
     - 300
     - T > 0
   * - Pressure
     - Pressure in bar to equilibrate during MM NPT step.
     - int
     - 1
     - P > 0
   * - MM minimization steps
     - Number of steps to minimize the MM potential energy of the system.
     - int
     - 2000
     - n > 0
   * - MM heat up steps
     - Number of steps to gradually increase the system's temperature with Langevin dynamics, allowing the simulation to reach a target temperature. Time step: 2 fs/step. Langevin dynamics collision constant: gamma_ln=2.0. Setting it to 0 skips the heating step. 
     - int
     - 10000
     - n >= 0
   * - MM NPT pressure equilibration steps
     - Number of steps to adjust the volume of the simulation box to reach a target constant pressure. Time step: 2 fs/step. Setting it to 0 skips the NPT step.
     - int
     - 300000
     - n >= 0
   * - MM NVE production run steps
     - Number of steps to evolve the system under constant particle number (N), volume (V), and energy (E) freely without any external constraints. Time step: 2 fs/step. Setting it to 0 skips the NVE step. 
     - int
     - 0
     - n >= 0
   * - QM/MM minimization steps
     - Number of steps to minimize the QM/MM potential energy of the system.
     - int
     - 250
     - n > 0
   * - QM/MM heat up steps
     - Number of steps to gradually increase the system's temperature with Langevin dynamics, allowing the simulation to reach a target temperature. Time step: 0.5 fs/step. Langevin dynamics collision constant: gamma_ln=5.0. Setting it to 0 skips the heating step. 
     - int
     - 250
     - n >= 0
  * - QM/MM NVE production run steps
     - Number of steps to evolve the system under constant particle number (N), volume (V), and energy (E) freely without any external constraints. Time step: 0.5 fs/step. Setting it to 0 skips the NVE step. 
     - int
     - 0
     - n >= 0
   * - QM method
     - Treating with high-level quantum mechanical accuracy. Available method: 'b3lyp', 'hf', 'case', 'dftb'.
     - str
     - b3lyp
     - NA
   * - Start Frame
     - First frame at which to start extracting from the trajectory the microsolvated clusters.
     - int
     - 0
     - n > 0
   * - Interval
     - Interval in frames at which to extract microsolvated clusters from the trajectory.
     - int
     - 100
     - n > 0
   * - Sell thickness
     - Solvent shell size for microsolvated clusters in Angstrom, upper limit for minimum solute-solvent distance.
     - int
     - 4
     - n > 0
