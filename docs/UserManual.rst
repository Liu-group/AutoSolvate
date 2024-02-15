Web-AutoSolvate Sudo-Parameter Definition
=========================================

Available job parameters are listed in the following Table.

.. list-table:: **General Parameter**
   :widths: 25 50 10 25 25
   :header-rows: 1

   * - Parameter
     - Description
     - Type
     - Default Value
     - Valid Range
   * - Solute
     - The molecule or group of molecules of interest. The file path to solute xyz file is accepted.
     - str
     - ""
     - Not given
   * - Solvent
     - The substance used to mimic the environment surrounding a solute molecule. Available solvents:'water', 'methanol', 'chloroform', 'nma', 'acetonitrile'.
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
     - Number of MM minimization steps.
     - int
     - 2000
     - n > 0
   * - MM heat up steps
     - Number of steps to increase the kinetic energy of the system, allowing the simulation to reach an equilibrium state and sample different conformations. Setting to 0 skips the MM heating step.
     - int
     - 10000
     - n > 0
   * - MM NPT pressure equilibration steps
     - Number of MM NPT steps, setting to 0 skips the MM NPT step.
     - int
     - 300000
     - n > 0
   * - MM NVE production run steps
     - Number of MM NVE steps, setting to 0 skips the MM NVE step.
     - int
     - 0
     - n > 0
