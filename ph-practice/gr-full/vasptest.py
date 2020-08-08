from pymatgen.io.vasp.inputs import Kpoints

k =  Kpoints.from_file('./KPOINTS')
print(k.kpts[0])