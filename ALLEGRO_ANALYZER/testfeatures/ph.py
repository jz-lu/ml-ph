# Testing use of phonopy API
from phonopy import Phonopy
from phonopy.structure.atoms import PhonopyAtoms
from phonopy.interface.calculator import read_crystal_structure
SCDIM = (6, 6, 1)
SCMATRIX = [[SCDIM[0], 0, 0], [0, SCDIM[1], 0], [0, 0, SCDIM[2]]]

# Build phonopy object from VASP unit cell POSCAR and generate displacements in supercell
uc, _ = read_crystal_structure("POSCAR-unit")
ph = Phonopy(uc, SCMATRIX, factor='VaspToTHz')
ph.generate_displacements()
scd = ph.supercells_with_displacements