import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import log
from pymatgen.io.vasp.inputs import Poscar # pylint: disable=import-error

a = [i for i in range(100)]
b = [i**2 for i in range(100)]
c = [log(i+1, 2) for i in range(100)]

p = Poscar.from_file("POSCAR")
lattices = p.as_dict()['structure']['lattice']['matrix'][:-1]
print(lattices)

x = [0, lattices[0][0], lattices[0][0] + lattices[1][0], lattices[1][0]]
y = [0, lattices[0][1], lattices[0][1] + lattices[1][1], lattices[1][1]]
z = np.array([i*i+j*j for j in y for i in x])
print(x)
print(y)
print(z)

X, Y = np.meshgrid(x, y)
Z = z.reshape(4, 4)

plt.pcolor(X, Y, Z, shading='auto')
# plt.fill(x, y, c="C0")
plt.savefig("test.png")

arr = np.array([a, b, c])
arr = np.transpose(arr)

