# Smoothen a noisy curve test
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter

def gaus(x, A=1, mu=0, sigma=1):
    return A * np.exp(-(x - mu)**2 / (2 * sigma**2))

A = 1; n = 101
x = np.linspace(-5, 5, n)
y = gaus(x, A=A)
jiggle = 0.1*(np.random.uniform(size=n)-0.5)
yerr = y + jiggle
ysmooth = savgol_filter(yerr, 31, 7) # window size 21, polynomial order 7
plt.plot(x, y, c='black', label="Pristine", alpha=0.9)
plt.plot(x, yerr, c='tab:green', label="Noisy", alpha=0.5)
plt.plot(x, ysmooth, c='red', label="Filtered")

plt.legend()
plt.show()

