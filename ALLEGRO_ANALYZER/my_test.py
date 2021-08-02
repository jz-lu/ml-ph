import phonopy
import numpy as np
import numpy.linalg as LA
from math import sqrt
from pymatgen.io.vasp.inputs import Poscar
import os

MAIN_DIR = '/Users/jonathanlu/Documents/' # must end with a '/'

def dsum(D, M):
    n = len(M)
    for i in range(n):
        s = np.split(D[3*i:3*(i+1)], n, axis=1)
        s = np.real_if_close(sum([sqrt(M[i]*M[j]) * blk for j, blk in enumerate(s)]))
        print(f"[At {i}] forces:\n{s}\n")

def fsum(f):
    sums = np.sum(f, axis=0)
    print(f"Sums:\n{sums}")
    
def adjust_force(f):
    print(f"FC TENSOR SHAPE: {f.shape}")
    sums = np.sum(f, axis=1)
    for i, s in enumerate(sums):
        f[i,i] -= s
        print(f"[{i}] Total change: {LA.norm(s)}")
    return f

p = Poscar.from_file(MAIN_DIR+"tmos2/POSCAR_LAYER1")
M = np.array(list(map(lambda x: x.atomic_mass, p.structure.species)))
assert M.shape == (3,)

os.chdir(MAIN_DIR+'tmos2/layer_1/analyses/phonon')
ph = phonopy.load(
    supercell_filename='SPOSCAR', 
    force_sets_filename="FORCE_SETS", 
    phonopy_yaml="phonopy_disp.yaml", 
    symprec=1e-4
    )
ph.produce_force_constants()
f = ph.force_constants
# fsum(ph.force_constants)
D = ph.get_dynamical_matrix_at_q([0,0,0])
print("Before:")
dsum(D, M)
print("After:")
fprime = adjust_force(f)
ph.force_constants = fprime
fpp = ph.force_constants
ph.set_force_constants(fprime)
assert np.allclose(fprime, fpp)
D = ph.get_dynamical_matrix_at_q([0,0,0])
dsum(D, M)

# print(f"D:\n{np.real_if_close(D)}")



