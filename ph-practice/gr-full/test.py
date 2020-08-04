from ___constants_names import *
from pymatgen.io.vasp.inputs import VaspInput, Incar, Poscar, Kpoints, Potcar
from pymatgen.io.vasp.outputs import Chgcar

chgcar = Chgcar.from_file('CHGCAR')
incar = Incar.from_file('INCAR')
poscar = Poscar.from_file('POSCAR')
potcar = Potcar.from_file('POTCAR')
kpoints = Kpoints.from_file('KPOINTS')
chgcar.write_file('TEST')

print(chgcar.data)
print(incar)

opDic = {CHGCAR_NAME: chgcar}
# vaspObj = VaspInput.from_directory(ROOT, opDic)
vaspObj = VaspInput(incar, kpoints, poscar, potcar, opDic)
print(vaspObj['CHGCAR'])

# vaspObj.write_input('./bobross')
# print(chgcar.as_dict())