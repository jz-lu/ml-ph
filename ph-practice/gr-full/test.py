from ___constants_vasp import *
from pymatgen.io.vasp.inputs import Kpoints

kpoints = Kpoints.from_file('../../pmg-practice/gr_BAND/KPOINTS')
print(kpoints)