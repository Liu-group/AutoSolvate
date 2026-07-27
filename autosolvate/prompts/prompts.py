prompt_system_type_choice = """AutoSolvate supports two system types based on your model: **Solute + Solvent** (Type 1) for a specific molecule(s) in a bulk environment, or **Homogeneous mixture** (Type 2) where all components are peers.

**1**. if you have one or a few special molecules (solutes) embedded in many other molecules (solvents).
**2**. if all components are present in comparable amounts, with no obvious "solute" to single out.

### examples (Type 1):
- *One citric acid molecule in water* (study hydration).
- *One perylene molecule in a methanol/cyclohexane mixture* (solute is perylene; solvents are the mixture).
- *One transition metal complex + counter ions in acetonitrile* (e.g., redox potential workflow).
- *A donor–acceptor pair in acetonitrile* (charge-transfer study).
- *A small cluster like 10 perylene molecules in a solvent mixture* (aggregation), where perylene is still the “solute set”.

### examples (Type 2):
- *An equimolar ethanol/water mixture* (miscibility, hydrogen-bond network).
- *A 30/70 acetonitrile/water mixture* (mixture dielectric/transport behavior).
- *A disordered box of tetracene molecules* (organic semiconductor morphology).
- *A mixed donor/acceptor material at similar fractions* (bulk electronic structure in disordered blends).
- *Electrolyte mixtures with multiple ionic species and solvent molecules* where no single “solute” is special.

### A quick decision rule
If you can point to a specific molecule (or small set) and say “**this is what I’m studying**,” choose **1**.
If you instead say “**the mixture itself is the system**,” choose **2**.
Always only type **1** or **2** to select your system type.

Please select the system type by typing **1** or **2** below:
> """

prompt_composition_method = """
Composition method:
1. Number. Specify the exact number of molecules for each component.
2. Weight portion. Specify the weight fraction for each component. 
3. Volume portion. Specify the volume fraction for each component. 
4. Molar portion. Specify the molar fraction for each component. 

Please select the composition method by typing **1**, **2**, **3**, or **4** below:
> """

prompt_single_component_method = """
Specify solvent amount by:
1. By molecule count (enter an exact integer).
2. By target density (g/cm^3) and molecular weight (g/mol); Then estimate molecule count from cube size.

Type **1** or **2** to choose.
> """

prompt_density_adjustment_choice = """
Density check options:
1. Suggest a new cube size to move the overall density toward the solvent density (preview only; you must confirm to apply).
2. Scale all solvent counts proportionally to the solvent density (preview only; you must confirm to apply).
3. Keep the current settings.

Type **1**, **2**, or **3** to choose.
> """

prompt_ask_solvent = """
How many solvent/components are present? (integer >=1 and <=10)
This refers to the number of distinct solvent types, not the total number of molecules.
> """

prompt_packmol_closeness = """
Packmol closeness controls minimum intermolecular distance.
- Default is 2.0 Å for general cases and custom solvents.
- For a single preset Amber solvent, we use its recommended closeness if available.

Press Enter to accept the suggested value, or type a new value in angstroms (e.g., 1.5):
> """

prompt_solvent_choice = """Select solvent/component {i}.

Eligible preset solvents (case-insensitive):
- water
- methanol
- chloroform
- nma
- acetonitrile

Type one preset name from the list above.

Type "custom" if:
- your solvent is not listed, OR
- you want to provide a custom structure file (.xyz/.pdb) and optionally your own frcmod/mol2/off/lib/prep.

Input examples:
- acetonitrile
- water
- custom
> """

prompt_n_solute_species = """How many distinct solute species? (integer >= 1)

Definition: a “solute species” corresponds to ONE structure file (.xyz/.pdb) that AutoSolvate will treat as a unit.

Examples:
- If you have two separate solutes A and B and you do NOT care about their initial relative geometry: provide 2 files -> enter 2.
- If your solute file contains multiple fragments / ion pairs that must keep a specific relative configuration: provide 1 multi-fragment file -> enter 1, then classify it as “multi-fragment/ion-pair”.

Type an integer (e.g., 1, 2, 3):
> """

prompt_classify_solute = """Classify this solute (choose 1/2/3):
1) Regular molecule (typical organic/ionic molecule)
2) Transition metal complex (requires MCPB + QM settings)
3) Multi-fragment / ion-pair complex (one file containing multiple fragments)

Tip: press Enter to accept the auto-detected suggestion shown in parentheses.
> """

prompt_ligand_chargefile = """Optional: ligand charge file path for MCPB (Amber MCPB charge file).

- Press Enter / type "skip" if you do not have one; AutoSolvate will try to determine ligand charges automatically.
- Otherwise, provide a file path.

> """

prompt_qm_software_choice = """Select QM software for MCPB.

Auto-detection will be printed above if available.
Type one of: orca / gau / g09 / g03 / gms
Press Enter to use the default (orca).

> """

prompt_qm_executable_path = """Path to QM executable.

- Press Enter / type "skip" to rely on your PATH (recommended if QM software is already on PATH).
- Otherwise, provide the full path to the executable (e.g., /path/to/orca).

> """

prompt_amberhome_path = """Provide AMBERHOME path.

- Press Enter / type "skip" to use the detected AMBERHOME shown above.
- Otherwise, provide the AMBERHOME directory path (e.g., /opt/amber22).

> """


prompt_final_confirmation_agent = """After reviewing the system setup above, 

Type **proceed** to start building the solvated system.
Type **write** to save the current configuration to a JSON file for later use and do not proceed.
Type **exit** to quit without saving or proceeding.
> """
# after determines copies
# check the solute type 
# then display the fragment
# display the number of copies for each residue:
# such as SUF: 1 copies, 5 atoms

# total charge

# the estimation is 1.0 g/cm3. the understanding is incorrect for some solvents such as chloroform

# if a preset solvent is used, and there is only a single solvent, set the default packmol closeness to that solvent's value