# ml-ph
Repository to calculate phonons for multi-layered two-dimensional materials

## Instructions on calling the automated script
You can run it directly or use a batch file to submit a job if you wish. The command-line calling, either way, should be of the form
```
python3 start.py <typeFlag: 0 for real space, 1 for config sampling> <I/O DIRECTORY> <vdW T/F> <GAMMA or MP> <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband, energies).\n\n\t If using input file...\t Usage: "python3 start.py -f <filename>" where parameters are in the input file separated by a newline
```
where _I/O DIRECTORY_ is the directory of the POSCAR and where the calculation will be run, _vdW (T/F)_ is a True/False on whether to use van der Waals interactions in the calculations or not, *GAMMA or MP* specifies the type of sampling on the Brillouin zone (for gamma-centered enter GAMMA or for Monkhorst-Pack enter MP), and the _calculation i_ specifies the types of calculations for electronic and phononic bands and DOS (inputs eledos, eleband, phdos, phband). Note that for _GAMMA or MP_, if you already have a valid KPOINTS file in your I/O directory, it does not matter which you put in; the file's sampling type will be read and used.

## Documentation
The general flow chart is given as a PDF in the docfiles folder, for configuration space sampling (the algorithm is similar for real-space calculations, omitting multiprocessing over the sampling and output related to displacement samples). The program can do both real space calculations and configuration space sampling based on type flag.

### Configuration Sampling Algorithm
Configurations are built as follows. Along with the standard INCAR, POTCAR, KPOINTS inputs (none of which are required—`Allegro Analyzer` automates the default construction of them all), one POSCAR must be given for each layer, for a minimum of two layers. These POSCARs are imported into _pymatgen_ objects in `Configuration::import_init_poscars()` in `__class_Configuration.py`. The lattice basis vectors, concatenated into a matrix, are extracted from the _pymatgen_ POSCAR object in `Configuration::get_lattices()`. A multilayered solid under the configuration sampling paradigm can differ _at most_ by a small constant scaling in the lattice; the lattice vectors themselves cannot be different up to normalization. Thus `Configuration::check_lattice_consistency()` subsequently examines the lattices over each layer and ensures that they are the same when normalized. If so, the calculation proceeds.


## Notes
Note that this code, at least the automated phonopy analysis part, likely does not work for 3-d materials or ones with  complicated band paths. For bulk materials the code can be modified easily on the phonopy analysis to use a package like seeKpath to automatically generate the high symmetry lines.

Any "special" notes:
* Input POSCARs should never begin with 'POSCAR-'. i.e. 'POSCAR_(stuff)' is fine and 'POSCAR' is good but not 'POSCAR-stuff'. The automated analysis relies on this keyword to analyze phonon structures properly.
