from pymatgen.io.vasp.inputs import Poscar

p = Poscar.from_file('POSCAR')
atoms_z_pos = []
atoms = p.structure.frac_coords

for i in atoms:
    if i[2] > 0:
        atoms_z_pos.append(round(i[2], 6))
atoms_z_pos = list(dict.fromkeys(atoms_z_pos))

interlayer_spacing = min(atoms_z_pos)
return interlayer_spacing