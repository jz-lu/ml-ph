import numpy as np
import matplotlib.pyplot as plt
from math import sqrt, log

def gaus(x, width):
    sigma = width / (2 * sqrt(2 * log(2))) 
    return np.exp(-np.power(x, 2.) / (2 * np.power(sigma, 2.)))

x = np.linspace(-10,10,100)
plt.plot(x, gaus(x, 10))
plt.show()