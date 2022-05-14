# Allegro Analyzer
Please read the general documentation in the main directory of `ml-ph` before reading this.

`Allegro Analyzer` consists of three pieces. The first is moderately fast, and involves the optimization of the lattice constant and interlayer spacing. The second is the slow piece, wherein VASP is used to compute the forces in configuration space. The third is the fast (and final part of the workflow) piece, wherein the forces are stitched together into the dynamical matrix. We discuss them separately below. Follow the workflow in the main documentation and refer to each piece as you get to it.

Throughout the workflow in the three pieces, you will have the option of using van der Waals corrections (via the R2SCAN-rVV10 functional) by specifying a `-v` flag. If you choose not to, the low-density approximation (LDA) will be used instead. No matter which you choose, hwoever, it is crucial that the choice is consistent throughout the entire workflow, as the lattice constants, interlayer spacings, GSFE coefficients, etc. differ significantly between the two.

## Constants (first piece)
By now you should have `POSCAR_LAYER1` and `POSCAR_LAYER2` ready to go in your directory of choice. `ml-ph` assumes that the monolayers POSCARs have been ionically relaxed from wherever you obtained them, e.g. the Materials project. If not, use VASP to relax them first. For each POSCAR, prepare a file `lc.txt` in the same directory as the POSCAR and in there specify the starting lattice constant, ending constant, and number of constants to uniformly sample in between; separate them by newlines. Next, run `start.py` with `-h` as a flag and familiarize yourself with the options. Most of them can be ignored, but it is important to specify that the type of calculation is `lc`. 

Once the calculations complete, run `optim_lc.py`, again using `-h` if needed. This will output the optimal lattice constant; that is, the lattice constant in the range specified in `lc.txt` that minimizes the total energy. There are two things to note.
1. In phonon calculations it is important that every pristine bilayer in configuration space is very well relaxed ionically. Thus the lattice constant should be precise up to third order (three decimal places). If you do not have such a good guess of precision, run a coarse sampling in a second-order range first, then run a third-order fine sampling around the optimal constant computed to second-order.
2. For homostructures you obviously only need to run the lattice constant optimization on one layer. For heterostructures, you must run them on both. Note that if the lattice constants differ by a second-order amount, it is likely not a stable bilayer and will produce nonsensical phonon calculations. Ensure that you choose heterostructures with similar lattice constants. For example, we found that MoSe2 and WSe2 had constants of 3.307 and 3.305 Angstroms, respectively, under the R2SCAN functional, which is a sufficiently small difference. Set the lattice constant for both materials to be the average of the two optimal constants.

You will also need to optimize the interlayer spacing. Although this is not strictly necessary, it will save many hours of computation time in the second piece. The calculation is similar, again using `start.py`, but this time the calculation type is `z` instead of `lc` and both POSCAR files are needed at once.

## DFT calculations (second piece)
Run `start.py` with calculation type `twist`. This will create three separate folders, `config`, `layer_1`, and `layer_2`. The former does the pristine bilayer calculations in configuration space and the latter does the monolayer calculations. You will need to specify some parameters in the command line, notably the interlayer spacing obtained from the first piece above. There are also optional flags, such as vdW functionals, that you can view in more detail with `-h`. By default the sampling will be `9 x 9` (`-s high`). For a low-sampling test run, use `-s low`. If you want to go higher, change `GRID_SAMPLE_HIGH` in `___constants_config.py`. Turn on `-f` (makes ionic relaxation more accurate). `e` specifies the initial cutoff energy for ionic relaxation. It doesn't really matter what this is for LDA, but for R2SCAN it can be hard to converge the relaxation for `e` too small, so try `1e-3` to start and if there are issues converging (see below) go bigger to `1e-2`. The calculation name `-n` is purely for cosmetic purposes if you want to give your SLURM job a specific name on the computing cluster. Use `--super 6` for a more accurate monolayer calculation.

Occassionally, calculations may fail because the ionic relaxation does not converge. You must resubmit these jobs until they do. Run `req.py --cfg` to check for and automatically resubmit failed jobs. Calculations of the forces may also fail for any number of reasons; run `req.py` to check. If you wish to skip over certain configurations because they are still running or on queue, use the `-s` flag and specify a list of the configuration indices.

Once complete, `cd` into `config` and run `config_analyze.py`, using `-h`. Then, record the output file `gsfe_coef.txt`, after checking that the fitted plot looks good. These are the coefficients `c0` through `c5` used in `Allegro Relaxer`.

## Phonon calculations (third piece)
All of the phonon calculations from the forces are performed by `twist_ph.py`. We outline band structure, density of states, and real space wavefunction.

Unless you are running a large angle, use `-r` so that the continuum relaxation is called. Otherwise, the bands will be absurd.

### Band structure
This is automatic just by calling `twist_ph.py`. Use `-h` for the details on the inputs. Use `-c` to specify the cutoff frequency in `1/cm`. Bands higher than the cutoff frequency will not be displayed.

### Density of states
Specify the `--dos` parameter in `twist_ph.py` (we recommend something around 45). This may be something you want to run on the computing cluster, as it will take several hours and many cores.

### Real space wavefunction
Provide input `rs.txt` which specifies which modes to output (starting at 0, the lowest energy mode). Provide input `k.txt` which specifies which k-points to do so at; each k-point should be specified as 3 numbers, separated by a space, on a single line. Put each k-point in a new line. You may specify either in Cartesian coordinates, or in direct coordinates of the reciprocal lattice basis. For the latter, add the `--kdir` flag when running `twist_ph.py`. 

Run `twist_ph.py` and include `--rs`. By default the plots are on a `2 x 2` moire supercell, but you can change this (e.g. if you are studying the 3-periodic point K) by `--rssz`.

### Theta space
If you want to see how the phonons at the Gamma point vary with twist angle, use `--thspc` as a flag. You will need to provide an input file `theta.txt` that gives the starting angle, ending angle, and number of angles to sample in between, separated by newlines. You also need `k.txt`, which is the same as in the real space calculation. This will output large `.npy` files containing the modes at each angle. You can plot whatever parts of the data relevant to your study. For more details, study `class ThetaSpaceDM` in `__class_DM.py`.