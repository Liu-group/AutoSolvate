AutoSolvateWeb Job Parameters
=============================

Available job parameters for the web interface are listed in the following Table.

.. list-table:: **General Parameter**
   :widths: 10 10 10 10 10
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
     - Spin multiplicity of solute.
     - int
     - 1
     - S >= 1
   * - Solute Cube Size
     - Size of the solvent box for Molecular Mechanics (Angstrom).
     - int
     - 54
     - L >=20
   * - Dry Run
     - Only generate the commands to run MD programs without executing.
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
     - Number of MM steps to minimize the potential energy of the system. Applicable to QM/MM calculation.
     - int
     - 2000
     - n > 0
   * - MM heat up steps
     - Number of steps to increase the kinetic energy of the system, allowing the simulation to reach an equilibrium state and sample different conformations. Setting it to 0 skips the heating step. Applicable to QM/MM calculation.
     - int
     - 10000
     - n >= 0
   * - MM NPT pressure equilibration steps
     - Number of steps to adjust the volume of the simulation box to maintain constant pressure while allowing the system to reach thermodynamic equilibrium. Setting it to 0 skips the NPT step. Applicable to QM/MM calculation.
     - int
     - 300000
     - n >= 0
   * - MM NVE production run steps
     - Number of steps to evolve the system under constant particle number (N), volume (V), and energy (E) freely without any external constraints. Setting it to 0 skips the NVE step. Applicable to QM/MM calculation.
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
