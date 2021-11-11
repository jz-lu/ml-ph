"""
Script for analyzing output of elastic moduli calculations kickedd off in `param.py`.
Python equivalent of `fit_elastic_const.m`
"""
import numpy as np
import argparse
import os
from param import folder_names
from gen_poscar_elastic import pathify
import matplotlib.pyplot as plt
from math import sqrt

parser = argparse.ArgumentParser(description="Output analyzer for elastic moduli calculations")
parser.add_argument("-d", "--dir", type=str, help="main working directory", default='.')
parser.add_argument("lconst", type=float, help="lattice constant")
args = parser.parse_args()
a0 = args.lconst
A0 = sqrt(3)/2 * a0**2

d = pathify(args.dir)
os.chdir(d)
infiles = [name + '.npy' for name in folder_names]
assert set(infiles).issubset(set(os.listdir())), f"Input files not found, run `read_elastic.py` first!"

bound = 0.04
eta = np.linspace(-bound, bound, 11)
Es = [] # as in energy E's
for infile in infiles:
    Es.append(np.load(infile))
Es = np.array(Es) # each row is all the energies for the given calc type across eta

# Do a quadratic fit on the energies, since elastic energy is quadratic in perturbation
fits = [0]*len(folder_names)
c2 = [0]*len(folder_names) # quadratic coefficient in fit
for i, E in enumerate(Es):
    z = np.polyfit(eta, E, 2)
    c2[i] = z[0]
    fits[i] = np.poly1d(z)
print("Successfully implemented quadratic fit")

# Plot the energies and fit as a function of displacement
eta_smooth = np.linspace(-bound, bound, 101)
for i, (E, name) in enumerate(zip(Es, folder_names)):
    plt.scatter(eta, E)
    plt.plot(eta_smooth, fits[i](eta_smooth), label=name)
plt.legend()
plt.savefig(d + "elastic_plots.png")
print(f"Saved plots to {d}")

# Compute moduli
gamma11 = c2[0]/A0 * 2
gamma12 = c2[1]/A0 - gamma11
gamma66 = c2[2]/A0 / 2
print(f"Coeffs: {c2}")

# in eV per cell 
G = gamma66 * A0
K = (gamma11 + gamma12)/2 * A0
print(f'G = {round(G, 4)} eV/cell')
print(f'K = {round(K, 4)} eV/cell')

