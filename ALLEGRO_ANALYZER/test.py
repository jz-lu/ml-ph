import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import log
import math
from pymatgen.io.vasp.inputs import Poscar # pylint: disable=import-error
import sys
from random import uniform as unif

p = Poscar.from_file("POS_TEST")
basis = p.as_dict()['structure']['lattice']['matrix'][:-1]
cob = np.transpose([i[:-1] for i in basis])
print("COB:", cob)

bset = np.arange(0, 3, 1)
pts = []

for i in bset:
    for j in bset:
        b = np.array([i, j])
        print("DIRECT:", b)
        c = list(np.dot(cob, b))
        c.append(i+j)
        print("CARTESIAN:", c)
        pts.append(c)

X = np.array([i[0] for i in pts])
Y = np.array([i[1] for i in pts])
z = np.array([i[2] for i in pts])
Z = z.reshape(int(math.sqrt(len(pts))), int(math.sqrt(len(pts))))
z2 = np.array([unif(0, 1) for i in pts])

print(pts)
print("X:", X)
print("Y:", Y)
print("z:", z)
print("z2:", z2)

fig, ax = plt.subplots()
cf = ax.tricontourf(X, Y, z2, levels=21, cmap="twilight_shifted")
fig.colorbar(cf, ax=ax)
ax.set_xlabel(r"$b_x$")
ax.set_ylabel(r"$b_y$")
ax.set_title(r"$E_{tot}(b)$")
fig.savefig("test.png")

