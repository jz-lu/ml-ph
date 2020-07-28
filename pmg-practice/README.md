# Notes on practice on graphene or bilayer

*  The lattice constant in this case is 2.4591. Note that the z-axis basis vector is not properly scaled to this but that's alright since it just needs to be large enough.
* In each folder, there is a Jupyter Notebook analysis.ipynb that conducts analysis.
* no_relax files use basic parsinng of pmg to get an energy vs spacing (interlayer for bilayer) plot.
* relax files are a relaxation with basic parsing.
* DOS and band structures analysis folders are labeled respectively. DOS uses same vasp calcualtions as relaxations, and  band structures has a different KPOINTS file (along high symmetry lines instead of mesh sampling).
