import numpy as np
from numpy.linalg import *
from phonopy import Phonopy
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
import os
import pymatgen.core.structure as struct
import pymatgen.io.vasp.inputs as inputs
import pymatgen.io.vasp.outputs as outputs
from math import *
import matplotlib.pyplot as plt
import sys
    
main = "/Users/jonathanlu/Documents/tmos2/layer_1/analyses/phonon/"
os.chdir(main)
ph = phonopy.load(unitcell_filename="POSCAR_unit", 
                   supercell_matrix=np.array([3,3,1]), 
                   primitive_matrix=np.eye(3),
                   log_level=1)
ph._log_level = 0

# unit cell 
mos2 = struct.Structure.from_file("POSCAR_unit")
A0 = mos2.lattice.matrix[0:2,0:2].T
G0 = 2*pi*inv(A0).T

# N x N x 1 supercell structure 
mos2_sc = struct.Structure.from_file("SPOSCAR")

nq = 101
Gamma = np.zeros(np.shape(G0[:,0]))
K = 2./3 * G0[:,0] + 1./3 * G0[:,1]
M = 0.5 * (G0[:,0] + G0[:,1])
corners = np.zeros([4,2])

corners[0, :] = K
corners[1, :] = Gamma
corners[2, :] = M
corners[3, :] = K
xticklabels = ('K', r'$\Gamma$', 'M', 'K')

k_set = np.zeros([nq*3-2,2]) # 3(nq-1) + 1 with kline[0] = kline[-1] for boundary 
k_mags = np.zeros([nq*3-2]) # flattening for plotting
kmag0 = 0
for kidx in range(np.shape(corners)[0]-1):
    idx_rng = kidx*(nq-1)
    kline_tmp = np.linspace(corners[kidx], corners[kidx+1], nq)
    if kidx != np.shape(corners)[0]-2:
        k_set[idx_rng:idx_rng+nq-1] = kline_tmp[0:nq-1]
    else:
        k_set[idx_rng:idx_rng+nq] = kline_tmp
    dk = corners[kidx+1,:] - corners[kidx,:]
    dk_norm = norm(dk)
    kmag_tmp = np.linspace(kmag0, kmag0+dk_norm, nq)
    kmag0 = kmag_tmp[-1]
    if kidx != np.shape(corners)[0]-2:
        k_mags[idx_rng:idx_rng+nq-1] = kmag_tmp[0:nq-1]
    else:
        k_mags[idx_rng:idx_rng+nq] = kmag_tmp
np.save("/Users/jonathanlu/Documents/k1.npy", k_mags)

corner_kmags = k_mags[0]
for kidx in range(np.shape(corners)[0]-1):
    corner_kmags=np.append([corner_kmags], k_mags[(kidx+1)*(nq-1)])

d_matrix = np.zeros([9, 9, k_set.shape[0]], dtype=complex)
evals = np.zeros([9, k_set.shape[0]], dtype=complex)
k_direct = k_set @ inv(G0).T
    
for q_idx, k in enumerate(k_set):
    d_matrix[:, :, q_idx] = ph.get_dynamical_matrix_at_q(k_direct[q_idx])
    evals[:, q_idx] = 15.633302**2 * np.real(np.sort(eigvals(d_matrix[:,:,q_idx])))

fig,ax=plt.subplots()
ax.plot(k_mags, np.transpose(np.sqrt(evals)), color='black')
ax.set_title("get_dynamical_matrix_at_q")
ax.set_xticks(corner_kmags)
ax.set_xticklabels(xticklabels)
ax.set_xlim([0, np.max(corner_kmags)])
plt.savefig("/Users/jonathanlu/Documents/ml-ph/ALLEGRO_ANALYZER/work.png")


