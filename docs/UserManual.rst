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
     - The molecule or group of molecules of interest, specifically structural and dynamic behavior over time in a solvent. The file path to solute xyz file is accepted.
     - str
     - ""
     - Not given
   * - Solvent
     - The substance used to mimic the environment surrounding a solute molecule in molecular interactions and dynamics studies. Available solvents:'water', 'methanol', 'chloroform', 'nma', 'acetonitrile'.
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
     - Pressure in bar to equilibrate during MM NPT step
     - int
     - 1
     - P > 0



`git <https://git-scm.com/>`_
