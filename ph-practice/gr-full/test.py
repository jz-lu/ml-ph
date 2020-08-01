from ___constants_vasp import *
from pymatgen.io.vasp.inputs import Kpoints

kpoints = Kpoints.gamma_automatic(kpts=(21, 21, 1), shift=(0,0,0))
kpoints.comment = 'jumpman jumpman always up to somethin'
print(kpoints)