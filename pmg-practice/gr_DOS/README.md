# Example of basic energy analysis, DOS on pmg

## VASP side
* First run a relaxation on the material and get a full calculation.
* Then set IBRION = -1 (forces only, no relaxation), ICHARG = 11 (keep charge densities constant the whole time to get consistent DOS), and NSW = 1, and run the VASP calculation again.
* Check analysis.ipynb for parsing of energies and DOS plotting example.