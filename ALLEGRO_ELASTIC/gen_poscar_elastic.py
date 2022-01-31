# -*- coding: utf-8 -*-
"""
Created on Mon Jul 15 12:50:58 2019

@author: zoe
Modified by: jonathan

Creating input file for elastic constant calculation
"""

import numpy as np
from math import sqrt
from pymatgen.core import structure
from pymatgen.io.vasp import inputs
from pymatgen.io.vasp.inputs import Poscar
import os, argparse, sys
import numpy.linalg as LA
import shutil

def pathify(path):
    path = os.path.abspath(path); assert os.path.isdir(path), f"Directory {path} does not exist"
    return path + ('' if path[-1] == '/' else '/')

# parser = argparse.ArgumentParser(description="Elastic moduli calculator for monolayers")
# parser.add_argument("mat", type=str, help="material to analyze, e.g. MoS2, Gr, etc.")
# # parser.add_argument("-d", "--dir", type=str, help="directory containing POSCAR", default='.')
# # parser.add_argument("-o", "--out", type=str, help="output directory", default='.')
# args = parser.parse_args()
# print(f"Args: {sys.argv[1:]}")

def get_poscar_elastic(args):
    ROOT = pathify(args.dir)
    mat = args.mat
    # fpath = ROOT + 'POSCAR'
    # poscar = Poscar.from_file(ROOT + 'POSCAR')
    # species= list(map(lambda x: x.species.elements[0].symbol, poscar.structure.sites))
    # lattice = poscar.structure.lattice.matrix
    # a0 = LA.norm(lattice[0])
    # langle = np.rad2deg(np.arccos(np.dot(lattice[0]/a0, lattice[1]/a0)))
    # assert np.isclose(langle, 60), f"Lattice angle must be 60, but is {langle}"

    print("WD:", ROOT)
    fpath = ROOT + mat + '/'
    if mat[0] == 'G':
        species = ['C', 'C']
        a0 = 2.457 # LDA: 2.446
        print(f"Material: Gr with lattice constants {a0}")
        pos = [[0, 0, 0], [sqrt(3)/3*a0, 0, 0]]
    elif mat[0] == 'W':
        species = ['W', 'Se', 'Se']
        a0 = 3.306 # LDA: 3.244
        print(f"Material:WSe2 with lattice constants {a0}")
        pos = [[0, 0, 0], [sqrt(3)/3*a0, 0, -1.671], [sqrt(3)/3*a0, 0, 1.671]]
    elif mat[:4] == 'MoSe':
        species = ['Mo', 'Se', 'Se']
        a0 = 3.306 # LDA: 3.246
        print(f"Material: MoSe2 with lattice constants {a0}")
        pos = [[0, 0, 0], [sqrt(3)/3*a0, 0, -1.565], [sqrt(3)/3*a0, 0, 1.565]]
    elif mat[0] == 'M':
        species = ['Mo', 'S', 'S']
        a0 = 3.178 # LDA: 3.122
        print(f"Material: MoS2 with lattice constants {a0}")
        pos = [[0, 0, 0], [sqrt(3)/3*a0, 0, -1.565], [sqrt(3)/3*a0, 0, 1.565]]
    elif mat[0] == 'J':
        species = ['Mo', 'S', 'Se']
        a0 = 3.1 # TODO
        print(f"Material: MoSSe (Janus) with lattice constants {a0}")
        pos = [[0, 0, 0], [sqrt(3)/3*a0, 0, -1.5569], [sqrt(3)/3*a0, 0, 1.5569]]

    if not os.path.isdir(fpath):
        os.mkdir(fpath)
    for car in ['INCAR', 'POTCAR', 'KPOINTS']:
        shutil.copyfile(ROOT+car, fpath+car)
    os.chdir(fpath)

    lattice = np.array([[sqrt(3)/2*a0, -0.5*a0, 0], [sqrt(3)/2*a0, 0.5*a0, 0], [0, 0, 20]])

    struct = structure.Structure(lattice, species, pos, coords_are_cartesian=True)
    struct.to(fmt='poscar', filename='POSCAR') # write unstreched POSCAR to file 
    N_ETA = 11

    eta = np.linspace(-0.04, 0.04, num=N_ETA)

    # stretching along the x direction 
    for i in range(len(eta)):
        try:
            os.mkdir('stretch1/')
        except: 
            pass
        os.chdir('stretch1/')
        try:  
            os.mkdir('def' + str(i) + '/')
        except: 
            pass 
        os.chdir('def' + str(i) + '/')
        lattice = np.array([[sqrt(3)/2*a0*(1+eta[i]), -0.5*a0, 0], [sqrt(3)/2*a0*(1+eta[i]), 0.5*a0, 0], [0, 0, 15]])
        struct = structure.Structure(lattice, species, pos, coords_are_cartesian=True)
        struct.to(fmt='poscar',filename='POSCAR')
        os.chdir('../../')
        
    # stretching along both x and y directions
    for i in range(len(eta)):
        try:
            os.mkdir('stretch2/')
        except: 
            pass
        os.chdir('stretch2/')
        try:  
            os.mkdir('def' + str(i) + '/')
        except: 
            pass 
        os.chdir('def' + str(i) + '/')
        lattice = np.array([[sqrt(3)/2*a0*(1+eta[i]), -0.5*a0*(1+eta[i]), 0], [sqrt(3)/2*a0*(1+eta[i]), 0.5*a0*(1+eta[i]), 0], [0, 0, 15]])
        struct = structure.Structure(lattice, species, pos, coords_are_cartesian=True)
        struct.to(fmt='poscar',filename='POSCAR')
        os.chdir('../../')  
        
    # shearing
    for i in range(len(eta)):
        try:
            os.mkdir('shear/')
        except: 
            pass
        os.chdir('shear/')
        try:  
            os.mkdir('def' + str(i) + '/')
        except: 
            pass 
        os.chdir('def' + str(i) + '/')
        lattice = np.array([[(sqrt(3) - eta[i])*a0*0.5, 0.5*a0*(-1+sqrt(3)*eta[i]), 0], [(sqrt(3) + eta[i])*a0*0.5, 0.5*a0*(1+sqrt(3)*eta[i]), 0], [0, 0, 15]])
        struct = structure.Structure(lattice, species, pos, coords_are_cartesian=True)
        struct.to(fmt='poscar',filename='POSCAR')
        os.chdir('../../')  
    return N_ETA, fpath
