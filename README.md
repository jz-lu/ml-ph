# ml-ph
This is `ml-ph` (`Allegro`), a scrappy library containing a set of tools that together compute the phononic physics of bilayered moire heterostructures with hexagonal pristine lattices. This library was developed for the purpose of a study at the Kaxiras Group at the Dept. of physics in Harvard University. Below provides some minimal documentation about the library, including the workflow for a phononic computation and the features available.

## How to cite this library
TODO

## Overciew
The calculation of phonons runs in two primary stages. The first involves computation of the pairwise forces of atoms in pristine bilayer configuration space over a discrete `N x N` set of configurations. (In the paper, we used `9 x 9`.) This step is built upon density functional theory (DFT) implemented by the Vienna _ab initio_ package (VASP), which requires separate licensing. See **Initial setup** for details. This calculation has been maximally parallelized but is still expected to take a long time. On average, we found this step to take about 8-24 hours per configuration to run on a computing cluster with 4000 GB RAM per core, 24 cores. DFT calculations will first require ionic relaxation, followed by force computations. The two steps are submitted as separate jobs. Be aware that since each configuration is relaxed, and then requires up to 36 force calculations in the frozen phonon method, the DFT calculations can total to a few thousand jobs submitted onto the computing cluster.

The second step involves the construction of the moire dynamical matrix from the forces. This step is very quick, and should take about 1-4 minutes for band structure and real space wavefunctions. However, for calculations of density of states (DOS), `ml-ph` must compute the dyamical matrix at every one of a `45 x 45` mesh in reciprocal space, which can take 2-10 hours. This step is based on the collective API calls of `phonopy`, `hiphive`, and `ml-ph` itself.

`ml-ph` has deprecated functions of calculating electronic band structures and DOS as well for pristine multilayered materials. We will not discuss how to do so there, but some of the old code is in `ALLEGRO_ANALYZER/__postprocess_relaxation.py` for the curious.

## Dependencies
Install Python packages via `conda` or `pip`. Look these up as needed for their own documentation on installation.

1. Python package `phonopy`, latest.
2. Python package `hiphive`, latest.
3. Vienna *ab initio* simulation package (VASP). Any version will do, as long as it is set up properly (see below).
4. A computing cluster that uses the SLURM scheduler. Most modern clusters use SLURM. For example, the Department of Energy NERSC computing clusters. (Without supercomputing power, this workflow is infeasible.)
5. Python 3.7.7 or later with the basics (`numpy`, `matplotlib`, etc.)
6. Python package `pymatgen`, latest.
8. Julia, any version at least `1.0`.

## Initial setup
The first order of business is to build VASP. This is standard to any condensed matter physics calculation and is not unique to `ml-ph`, so we will not discuss it in detail. Whatever the method is, update `ALLEGRO_ANALYZER/STATIC_BAT_DNE` with the right command to call the VASP executable. 

Ensure next that the python path and the location of the VASP POTCAR (pseudopotential files)are correctly set up in `~/.pmgrc.yaml`. The details of this can be found in the `pymatgen` setup documentation.

Depending on how different your implementation of SLURM is from ours, you may need to do significant modification to `ALLEGRO_ANALYZER/make_exescr.py`, which builds the general bash executable script that calls SLURM. See that file for details.

## Workflow
Every VASP calculation requires the INCAR, KPOINTS, POSCAR, and POTCAR input files. See VASP documentation about these. `ml-ph` is set up to automatically create the INCAR, KPOINTS, and POTCAR files for you, but you must provide a POSCAR file in the directory wherein you run the calculations. If you want to specify your own INCAR/KPOINTS, you may add that as an input file along the required POSCAR, and `ml-ph` will use that in lieu of a default creation. If you want to change the parameters of the default creation, check `ALLEGRO_ANALYZER/___constants_vasp.py`.

The general workflow is as follows. Prepare two POSCAR files, `POSCAR_LAYER1` and `POSCAR_LAYER2` for the two monolayers. Ensure that they are positioned such that if they were stacked on top of each other directly, they would be at the highest energy (AA) stacking configuration, for `ml-ph` assumes this. Moreover, the vacuum separation (the length of the out-of-plane lattice basis vector) is assumed to be 35 for `ml-ph` calculations. This can be changed in `ALLEGRO_ANALYZER/___constants_vasp.py`. However, it is important that the vacuum separation is sufficiently large, so it is better to change your POSCARs to match the assumed 35.0. `ALLEGRO_ANALYZER/modpos.py` does this automatically; run it with `-h` for details.

Each part of the workflow utilizes a different subroutine of `Allegro`, and the details are found in the `README.md` in the respective folders.

1. Obtain the elastic constants `K` and `G` of each monolayer by calling `ALLEGRO_ELASTIC`.
2. Run a twisted calculation in `ALLEGRO_ANALYZER`. Obtain the GSFE coefficients `c0` through `c5` as described in the documentation there. For phonon calculations, you may discard `c0`.
3. Create a parameters file in `ALLEGRO_RELAXER` for your material. Use the templates in the example files. Copy in the `K` and `G` for each layer as well as `c1` through `c5`.
4. Get the band structure, real space wavefunction, and density of states to your heart's content by running some of the scripts in `ALLEGRO_ANALYZER`.
5. Go forth and prosper!

See Fig. 1 in the paper from **How to cite** for a visualization of the workflow.