# Dynamical matrix calculations for twisted materials
import numpy as np
import sys, copy, os
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
import numpy.linalg as LA
from math import pi, floor
from itertools import product as prod
from ___helpers_parsing import succ, warn, err, is_flag, check_not_flag

USAGE_MSG = 'Usage: python3 <DIR>/dm.py -deg <twist degree> -dir <main dir> -o <output dir>'

def G_A_moire(theta, log=False):
    # Generate G basis via rotation matrices
    R1 = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]])
    R2 = LA.inv(R1)
    G1 = np.matmul(R1, G0)
    G2 = np.matmul(R2, G0) 
    GM = G1 - G2 # moire G-basis
    AM = LA.inv(GM).T / (2 * pi) # moire A-basis
    if log:
        print("Moire G-basis:\n", GM)
        print("Moire A-basis:\n", AM)
    return (GM, AM)

class BZSampler:
    def __init__(self):
        self.g_idxs = None; self.g_arr = None
        return
    # Sample Gtilde vectors
    def sample_G(self, G_moire, mn_grid_sz, max_shell):
        grid = np.arange(-mn_grid_sz, mn_grid_sz + 1)
        GM1 = G_moire[:,0]; GM2 = G_moire[:,1] # Moire reciprocal lattice vectors
        g_arr = np.array([np.array([m*GM1 + n*GM2]) for m, n in prod(grid, grid)])
        g_idxs = np.array(list(prod(grid, grid)))

        # Filter out any G that does not satisfy the closeness condition |GM| < (shell+0.1) * |GM1|
        g_cutoff = (floor(max_shell)+0.1) * LA.norm(GM1) # add 0.1 for numerical error
        cut_makers = (np.sqrt(g_arr[:,0]**2 + g_arr[:,1]**2) <= g_cutoff) # indicators for which G made the cutoff
        g_idxs = g_idxs[cut_makers,:]
        g_arr = g_arr[cut_makers,:] # special numpy syntax for filtering by an indicator array `cut_makers`
        self.g_idxs = g_idxs; self.g_arr = g_arr
        return (g_idxs, g_arr)

    # Sample k points along IBZ boundary line
    def sample_k(self, nk, G0):
        G00 = G0[:,0]; G01 = G0[:,1]; d = G00.shape
        Gamma = np.zeros(d); K = 1/3 * (G00 + G01); M = 1/2 * G00
        ptlen = 4; pt = np.zeros([ptlen, d]) # k-points boundaries
        pt[0,:] = Gamma; pt[1,:] = K; pt[2,:] = M; pt[3,:] = Gamma
        for kidx in range(ptlen-1):
            pass # TODO
        return
    
    def plot_sampling(self):
        assert self.g_idxs and self.g_arr
        # TODO
        return

# Compute dynamical matrix block element for given q = Gtilde and phonopy object ph
def dm(q, ph):
    pass #! import

# Create level-1 block matrix
def block_l1():
    pass # TODO

# Create level-2 block matrix with intralayer and interlayer terms
def block_l2(D_intras, D_inter):
    assert len(D_intras) == 2
    return np.block([[D_intras[0], D_inter], [D_inter.conjudate().T, D_intras[1]]])

# Parse input
args = copy.deepcopy(sys.argv)[1:]; i = 0; n = len(args)
theta = None; indir = '.'; outdir = None
while i < n:
    if not is_flag(args[i]):
        warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
        i += 1; continue
    if args[i] == '-deg':
        i += 1; check_not_flag(args[i]); theta = np.deg2rad(float(args[i])); i +=1
    elif args[i] == '-dir':
        i +=1; check_not_flag(args[i])
        indir = args[i]
        if args[i] in ['.', './', '..', '../'] or args[i][0] == '.':
            warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
        i +=1
    elif args[i] == '-o':
        i +=1; check_not_flag(args[i])
        outdir = args[i]
        if args[i] in ['.', './', '..', '../'] or args[i][0] == '.':
            warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
        i +=1

if not (theta and indir):
    err(USAGE_MSG)
elif 180 < theta < 0:
    err(f"Error: invalid twist angle {theta} degrees")
elif not outdir:
    outdir = indir
assert os.path.isdir(indir) and os.path.isdir(outdir)

# Build realspace A-basis and reciprocal space G-basis in the moire cell
s = struct.Structure.from_file("POSCAR_MoS2")
A0 = s.lattice.matrix[0:2, 0:2].T # makes realspace lattice A-basis matrix [a1 a2]
G0 = 2 * pi * LA.inv(A0).T 
print("A0:\n", A0)
print("G0:\n", G0)
