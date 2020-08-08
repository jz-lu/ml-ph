# ml-ph
Repository to calculate phonons for multi-layered two-dimensional materials

## Instructions on calling the automated script
You can run it directly or use a batch file to submit a job if you wish. The command-line calling, either way, should be of the form
```
python3 start.py <I/O DIRECTORY> <vdW (T/F)> <GAMMA or MP> <calculation 1> ... <calculation n>
```
where _<I/O DIRECTORY>_ is the directory of the POSCAR and where the calculation will be run, _<vdW (T/F)>_ is a True/False on whether to use van der Waals interactions in the calculations or not, _<GAMMA or MP>_ specifies the type of sampling on the Brillouin zone (for gamma-centered enter GAMMA or for Monkhorst-Pack enter MP), and the _<calculation i>_ specifies the types of calculations for electronic and phononic bands and DOS (inputs eledos, eleband, phdos, phband).


## Notes
Note that this code, at least the automated phonopy analysis part, likely does not work for 3-d materials or ones with  complicated band paths. For bulk materials the code can be modified easily on the phonopy analysis to use a package like seeKpath to automatically generate the high symmetry lines.

Any "special" notes:
* Input POSCARs should never begin with 'POSCAR-'. i.e. 'POSCAR_(stuff)' is fine and 'POSCAR' is good but not 'POSCAR-stuff'. The automated analysis relies on this keyword to analyze phonon structures properly.
