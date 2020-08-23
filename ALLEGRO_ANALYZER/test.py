import numpy as np
import scipy.linalg as la
import sys

from pymatgen.io.vasp.inputs import Poscar

p = Poscar.from_file('POSCAR')
sd = []
at = p.structure.pop(1)
at.frac_coords = at.frac_coords + np.array([5, 5, 1.3])
at.frac_coords = at.frac_coords % 1
print(at.species, at.frac_coords)
p.structure.append(at.species, at.frac_coords)
print(p.structure.num_sites)
for i in range(p.structure.num_sites):
    sd.append([False, False, True])
p = Poscar(p.structure, selective_dynamics=sd)
print(p)
sys.exit()

# p = Poscar.from_file('POSCAR')
# print(dir(p.structure))
# print(p.structure.frac_coords)
# mat = p.structure.frac_coords
# for i in range(len(mat)):
#     print(mat[i][2])

# p_dict = p.as_dict()
# lattice = np.array(p_dict['structure'])
# print(lattice)



# print(np.linalg.norm(lattice[0]), np.linalg.norm(lattice[1]), np.linalg.norm(lattice[2]))
# print(p_dict['structure'])
# print(np.where((lattice[:][2] == 0))) # Get which rows in lattice vector-columns matrix have zero z-coord
# print(lattice[2][2])

# print(np.linalg.norm(lattice[0]) - np.linalg.norm(lattice[1]))
# cob_matrix = np.transpose(lattice) # Make lattice vectors columns for change of basis

# print(cob_matrix)

# test = np.array([1, 0, 0]).reshape(3, 1)
# print(test)

# print('product:\n', cob_matrix @ test)