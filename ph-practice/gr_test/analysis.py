import phonopy 
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr_test'

# structure = Structure.from_file('POSCAR_old')
# print(structure)

# ibz = HighSymmKpath(structure)
# print(ibz)

# kpts = Kpoints.automatic_linemode(50, ibz)
# print(kpts)

# kfile = Kpoints.from_file(ROOT + '/LINE_KPOINTS')
# kfile.comment = 'C2-line'
# print(kfile)

poscar = Poscar.from_file('POSCAR')
print(poscar.natoms[0])