# ml-ph
Repository to calculate phonons for multi-layered two-dimensional materials

_Documentation Notes_:
1. This documentation regularly switches between the terminology of the "fixed" layer and the "first" layer (i.e. bottom layer). They both refer to the layer serving as the reference point for the configuration; this layer sets the standard of the lattice constant that all other layers will be strained to in the strain-shift algorithm. Since the choice of the fixed layer is arbitrary, we will always take it to be the first or bottommost layer.
2. "Selective dynamics", which refers to a constraint on which directions (in terms of the lattice basis vector set) the VASP DFT calculator's relaxation protocol is allowed to relax the solid on, is shortened to "SD".

## Overview
The purpose, (blackbox style) functions, and sample calculations are found in the _docs_ folder. The below details how to run the script, and then discusses a high-level idea of how it works. **A flow chart detailing the algorithm is found in the _docs_ folder as well**.


## Instructions on calling the automated script
If any of this section if unclear, read the documentation first, in particular the **Input** section. You can run it directly or use a batch file to submit a job if you wish. The command-line calling, either way, should be of the form
```
python3 start.py <typeFlag: 0 for real space, 1 for config sampling> <I/O DIRECTORY> <vdW T/F> <GAMMA or MP> <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband, energies).\n\n\t If using input file...\t Usage: "python3 start.py -f <filename>" where parameters are in the input file separated by a newline
```
where _I/O DIRECTORY_ is the directory of the POSCAR and where the calculation will be run, _vdW (T/F)_ is a True/False on whether to use van der Waals interactions in the calculations or not, *GAMMA or MP* specifies the type of sampling on the Brillouin zone (for gamma-centered enter GAMMA or for Monkhorst-Pack enter MP), and the _calculation i_ specifies the types of calculations for electronic and phononic bands and DOS (inputs eledos, eleband, phdos, phband). Note that for _GAMMA or MP_, if you already have a valid KPOINTS file in your I/O directory, it does not matter which you put in; the file's sampling type will be read and used.

## Documentation
The general flow chart is given as a PDF in the docfiles folder, for configuration space sampling (the algorithm is similar for real-space calculations, omitting multiprocessing over the sampling and output related to displacement samples). The program can do both real-space calculations and configuration-space sampling based on type flag. The following 

### Input
For any calculation, the four inputs are the INCAR (VASP settings), POTCAR (pseudopotential), KPOINTS (reciprocal space sampling line/grid), and POSCAR (solid spcification). If running a real-space calculations, the input is a single POSCAR; for configuration space, the input is one POSCAR per layer (see **Notes** on input format) and the fixed layer is the file that comes *first in alphanumeric order*. INCAR, KPOINTS (grid for DOS), and POTCAR are optional—`Allegro Analyzer` will generate them with default options if they are not found. For band structure, a line KPOINTS file named `LINE_KPOINTS` should be inputted. All inputs should go into the I/O directory specified on the command line.

### Configuration Sampling

#### Input validation
Configurations are built as follows. Along with the standard INCAR, POTCAR, KPOINTS inputs (none of which are required—`Allegro Analyzer` automates the default construction of them all), one POSCAR must be given for each layer, for a minimum of two layers. These POSCARs are imported into _pymatgen_ objects in `Configuration::import_init_poscars()` in `__class_Configuration.py`, from individual files that must be located at the specified `ROOT` I/O directory from the command line. The lattice basis vectors, concatenated into a matrix, are extracted from the _pymatgen_ POSCAR object in `Configuration::get_lattices()`. The fixed layer (which is the first layer under our model) in-plane lattice vectors (the 2 basis vectors in the plane of the 2D solid) are checked to have the same norm (they should be scaled identically), and then normalized in `Configuration::__get_normed_fixed_lattice()`. A multilayered solid under the configuration sampling paradigm can differ _at most_ by a small constant scaling in the lattice; the lattice vectors themselves _must be identical_ up to normalization. Thus `Configuration::check_lattice_consistency()` subsequently examines the lattices over each layer and ensures that they are the same when normalized. Finally,  If any of the above checks fail, `Allegro Analyzer` terminates. 

#### POSCAR construction and the strain-shift algorithm
Validated input POSCARs, one for each layer, are concatenated into a single POSCAR as follows. Shifts are generated via the static class method `Configuration::sample_grid()`. Given a shift, a single POSCAR containing all the layers are generated via `Configuration::build_config_poscar()`, which implements the strain-shift (SS) algorithm. Informally, SS strains every layer except the fixed so that every layer has the same lattice basis (we assume that this strain maintains solid stability, so that the lattice constant mismatch is sufficiently small, as explained in **Input validation**), then shifts every sublattice atom modulo the unit cell torus. SD is allowed only in the interlayer direction for all but the fixed layer. In practice this simply amounts to a few lines of code:
```
for i in nonfixed_layers:
  for j in i->atoms:
    j + shift_vec mod (1, 1, 1)
```
No explicit straining is necessary since the fixed layer serves as the base POSCAR, and each atom from the other layers are subsequently added to the base POSCAR after the shift, so by virtue of sharing the same lattice basis as the fixed, the strain has already been made.

### DFT Relaxation and Energy Computation
The program assumes that the individual layers are either already relaxed, or are strained anyway so it isn't actually supposed to be relaxed (e.g. mismatched lattices in multilayers, twists, etc.). Thus only the interlayer spacing is relaxed. The relaxation process is computed with DFT handled by the program VASP (Vienna Ab Initio Simulation Package). The inputs required on VASP are the [INCAR](https://www.vasp.at/wiki/index.php/INCAR) (VASP settings), [POSCAR](https://www.vasp.at/wiki/wiki/index.php/POSCAR) (gives all specs of the unit cell of the solid), [POTCAR](https://www.vasp.at/wiki/wiki/index.php/POTCAR) (which is always handled entirely by the program and does not require user knowledge--it is just a file giving the pseudopotential that VASP is to use in simluating the given material), and [KPOINTS](https://www.vasp.at/wiki/wiki/index.php/KPOINTS) (specifies the sampling of the BZ, a mesh/grid for DOS and a line around the IBZ for band structure). To relax the solid, VASP starts by assuming a trial wavefunction using Linear Combination of Orbitals (LCAO), then iteratively solves the corresponding Schrodinger equation until a self-consistent wavefunction is found (the calculation is said to *converge* when this occurs). For phonon calculations, the convergence threshold should be extremely small for accurate calculations, and in the default INCAR it is set to `10^-8` energy units. In each round, the energy computed, along with use of the Hellmann-Feynman theorem (get the forces from the Hamiltonian matrix elements), give the direction and magnitude by which VASP nnudges the atoms, thus "relaxing them" to an equilibrium point. Overall, the process is 2-layered: there are many rounds of relaxation where the atoms are moved based on the computed forces until an equilibrium threshold, set in INCAR, has been reached (or failure to relax is outputted), and in each relaxation round the wavefunctions are solved iteratively until self-consistent and the forces computed. The total and Fermi energy are then outputted by the program.

### Electronic Calculations
TODO

### Phononic Calculations
TODO

## Notes
Note that this code, at least the automated phonopy analysis part, likely does not work for 3-d materials or ones with  complicated band paths. For bulk materials the code can be modified easily on the phonopy analysis to use a package like seeKpath to automatically generate the high symmetry lines.

Any "special" notes:
* Input POSCARs should never begin with 'POSCAR-'. i.e. 'POSCAR_(stuff)' is fine and 'POSCAR' is good but not 'POSCAR-stuff'. The automated analysis relies on this keyword to analyze phonon structures properly.

## Informal Changelog (major only)
(Use `git log --oneline --decorate --color` for detailed changelog).
2021-01-17: fixed strain bugs in Configuration class and fixed z-lattice vector constant off by factor of 10.
