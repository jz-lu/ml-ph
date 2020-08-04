from ___constants_names import *
from pymatgen.io.vasp.inputs import VaspInput
from pymatgen.io.vasp.outputs import Chgcar

vaspObj = VaspInput.from_directory(ROOT, {CHGCAR_NAME: Chgcar})
vaspObj.write_input('./bobross')