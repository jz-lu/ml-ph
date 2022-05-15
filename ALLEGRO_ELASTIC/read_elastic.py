# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 18:11:46 2019

@author: zoe
Modified by: jonathan
"""
import numpy as np
import matplotlib.pyplot as plt
import os, argparse
from pymatgen.core import structure
from pymatgen.io.vasp import outputs
from param import folder_names
from gen_poscar_elastic import pathify

parser = argparse.ArgumentParser(description='Parse results of elastic moduli calculations')
parser.add_argument('-d', '--dir', type=str, help='working directory', default='.')
arg = parser.parse_args()
main_dir = pathify(arg.dir)

for folder in folder_names:
    fpath = main_dir + folder
    assert os.path.isdir(fpath)
    ind = 0; E0 = []; a1x = []; a1y = []; a2x = []; a2y = []; xf = []; yf = []
    while os.path.isfile(fpath + '/def' + str(ind) + '/POSCAR'):
        # read poscar
        struct = structure.Structure.from_file(fpath + '/def' + str(ind) + '/CONTCAR')
        pos = struct.cart_coords
        matrix = struct.lattice.matrix
        xf.append(pos[1][0])
        yf.append(pos[1][1])
        a1x.append(matrix[0][0])
        a1y.append(matrix[0][1])
        a2x.append(matrix[1][0])
        a2y.append(matrix[1][1])
        osz = outputs.Oszicar(fpath + '/def' + str(ind) + '/OSZICAR')
        E0.append(float(osz.final_energy))
        ind = ind + 1

    writepath = main_dir + folder + '.txt'
    apath = main_dir + folder + '.npy'
    with open(writepath, "a+") as f:
        f.write('a1x, a1y, a2x, a2y, xf, yf, E0\n')
        for i in range(len(xf)): 
            f.write(str(a1x[i]) + ', ' + str(a1y[i]) + ', ' + str(a2x[i]) + ', ' + str(a2y[i]) + ', ' + str(xf[i]) + ', ' + str(yf[i]) + ', ' + str(E0[i]) + '\n')
    np.save(apath, E0)
    print(f"Finished writing to {writepath}")
print("All analyses completed. Run `fit_elastic_const.py` to obtain moduli coefficients.")

