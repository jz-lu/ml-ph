# Allegro Elastic
Please read the general documentation in the main directory of `ml-ph` before reading this.

This is `Allegro Elastic`, an arm of `Allegro` responsible for the calculation of the elastic constants `K` (bulk modulus) and `G` (shear modulus) for monolayered materials.

Regrettably, `Allegro Elastic` is less well-designed than the other arms of `Allegro`. In particular, before you compute anything, you will have to add the material you want in the branches of `get_poscar_elastic()` in `gen_poscar_elastic.py`. Give it a name, and add the relevant information by following the template from materials we already have, including graphene, MoS2, WSe2, etc. If the monolayers of interest are already there, simply adjust the lattice constant `a0`. 

Place the `POSCAR` file in the directory where you wish to perform the calculation and call `param.py`. For the elastic calculation, you will need to give the INCAR, POTCAR, and KPOINTS files. Use a moderately large k-point sampling (e.g. `15 x 15`). The INCAR should be a non-relaxing one (i.e. electronic calculations only, such as the one produced in the DFT force calculation step from `Allegro Analyzer`). The POTCAR may be copied from the `Allegro Analyzer` calculations.

`param.py` will compute all the relevant shearing and stress constants in the computing cluster. You will need to adjust `param.py` first, by changing the form of the SLURM submission file to match that of your computing cluster. 

Once the computing cluster is complete, `cd` into the directory containing all of the calculations and run `read_elastic.py`. Subsequently, run `fit_elastic_const.py`. All of these have the `-h` function on the command line, which prints their documentation and required input paramters. This will output `K` and `G` as well as a plot of the fit done to create it. If it looks smooth and parabolic, the calculation went well and you can trust the elastic constants.