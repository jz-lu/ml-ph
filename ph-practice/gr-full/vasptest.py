from pymatgen.io.vasp.inputs import VaspInput

vasp_obj = VaspInput.from_directory('.')
print(vasp_obj['INCAR'])
vasp_obj.run_vasp()