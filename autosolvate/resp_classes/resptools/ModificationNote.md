2023-09-24

The file, 

python3_openbabel/resp_gamess.py, 

was modified from resptools/python/Resp.py of the resptools repository 
https://github.com/choderalab/resptools developed by John Chodera 

Changes were made to make the module fit the use cases in AutoSolvate.

1. Use python3 syntax
2. Removed all dependencies on OpenEye, but use OpenBabel
2. Added support for open-shell molecules (spin multiplicity>1)
3. Only consider single-configuration RESP fitting. Remove all compoents related to multi-component resp fitting
4. Job control with SLURMS instead of SGE/PBS
5. Changes to make the processing of GAMESS output file match GAMESS VERSION = 31 JUL 2022 (R1)
