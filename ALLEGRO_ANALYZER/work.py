import numpy as np
import matplotlib.pyplot as plt
from itertools import product as prod

A = np.array([[3.12, 0], [.561000, 2.703732]]).T
line = np.linspace(0, 1, 11)
assert line.shape == (11,)
grid = np.array(list(prod(line, line))) @ A.T
x = grid[:,0]; y = grid[:,1]
Z = np.random.rand(11*11)
# x = np.arange(cols)
# y = np.arange(rows)
plt.tricontourf(x, y, Z, levels=301)
plt.colorbar()

plt.show()