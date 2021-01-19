import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from math import log

a = [i for i in range(100)]
b = [i**2 for i in range(100)]
c = [log(i+1, 2) for i in range(100)]

arr = np.array([a, b, c])
arr = np.transpose(arr)

with open("foo.csv", 'wb') as f:
    f.write(b"a, b, c\n")
    np.savetxt("foo.csv", arr, delimiter=",")
