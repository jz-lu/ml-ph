import phonopy 
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from pymatgen.core.structure import Structure
from pymatgen.symmetry.bandstructure import HighSymmKpath
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr_test'

k = Kpoints.from_file('LINE_KPOINTS')
print(k.kpts)