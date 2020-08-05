from pymatgen.io.vasp.inputs import VaspInput, Incar, Kpoints

vasp_obj = VaspInput.from_directory('.')

potcar = vasp_obj['POTCAR']
print(potcar.symbols)
print(vasp_obj['KPOINTS'])