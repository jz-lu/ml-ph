# ml-ph
Repository to calculate phonons for multi-layered two-dimensional materials

Note that this code, at least the automated phonopy analysis part, likely does not work for 3-d materials or ones with  complicated band paths. For bulk materials the code can be modified easily on the phonopy analysis to use a package like seeKpath to automatically generate the high symmetry lines.

Any "special" notes:
* Input POSCARs should never begin with 'POSCAR-'. i.e. 'POSCAR_(stuff)' is fine and 'POSCAR' is good but not 'POSCAR-stuff'. The automated analysis relies on this keyword to analyze phonon structures properly.
