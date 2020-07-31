# Guide to generating a phonon DOS

* Before starting, ensure phonopy is installed properly and that you create an environment for it, and that the ennvironment is activated. Details about it can be found in the section *Minimum steps to install and use phonopy via conda* in [this link](https://phonopy.github.io/phonopy/install.html).

**DOS**
* (Assuming you have used a grid mesh for KPOINTS and thus are doing a DOS calculation, see the second half for band calculations) Use the standard input files to start out.
* Compute the relaxation for the structure with a fine EDIFF = 1E-8 and EDIFFG the same. Such levels of accuracy are necessary avoid negative frequencies of phonon bands.
* Once geometric relaxation is complete, run `phonopy -d --dim="X Y Z" -c POSCAR-unit`, where X Y Z are the supercell sizes in each direction of the PUC vectors and POSCAR-unit is the name of the unit cell POSCAR. This will create a number of POSCAR files POSCAR-XXX, where XXX is some number iterating from 001, depending on the symmetry of the material. They represent the atomic displacements.
* Run a VASP force calculation (IBRION = -1) for each POSCAR-XYZ, in separate subfolders. ISIF = 2 (or remove it, it is 2 by default).
* Run `phonopy -f disp-001/vasprun.xml disp-002/vasprun.xml disp-003/vasprun.xml` and so on assuming you have named the subfolders disp-XYZ accordingly. This creates FORCE_SETS for each displacement (one FORCE_SETS for all displacements, so they have to be run together). For more concise notation, run `phonopy -f disp-{001..XYZ}/vasprun.xml` instead.
* Create a new file `mesh.conf` in the original folder with the following specifications: ```
ATOM_NAME = LIST ATOMS IN PUC SEPARATED BY A SPACE
DIM = SPECIFY SUPERCELL DIM SAME AS PREPROCESSING --dim="X Y Z" FLAG
MP = SPECIFY MESH (IF UNSURE JUST PICK SAME AS KPOINTS MESH)
```
* Run `phonopy -p mesh.conf -c POSCAR-unit` to get the DOS. Do this on an interactive system (local or interactive session) to get an output image file, or else it will output a simple text-style file that is easy to parse and plot.

**Band Structure**
* Do the same thing as before, setting KPOINTS to a line appropriate to the reciprocal lattice IBZ boundary.
* Instead of creating `mesh.conf`, create `band.conf` with the following specifications:
```
ATOM_NAME = SAME AS BEFORE
DIM =  SAME AS BEFORE
BAND = LINE SPECIFICATIONS SAME AS IN KPOINTS, e.g. 0.5 0.5 0.5  0.0 0.0 0.0  0.5 0.5 0.0  0.0 0.5 0.0
```
* Run `phonopy -p band.conf -c POSCAR-unit`.
