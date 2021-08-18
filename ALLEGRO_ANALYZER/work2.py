import numpy as np
from numpy.linalg import *
from phonopy import Phonopy
import phonopy
from phonopy.phonon.band_structure import get_band_qpoints_and_path_connections
import os
import pymatgen.core.structure as struct
import pymatgen.io.vasp.inputs as inputs
import pymatgen.io.vasp.outputs as outputs
import h5py
from math import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import sys

def sc_idx(sc, uc):
    uc_coords = uc.cart_coords; sc_coords = sc.cart_coords
    sc_nat = len(sc_coords); uc_nat = len(uc_coords) # num sc/uc atoms
    pos_sc_id = []; n_uc = int(sc_nat/uc_nat) # num unit cells

    # SPOSCAR arranges all atoms of each type contiguously, so the indexes must
    # be the same for each contiguous region of `n_uc`.
    for i in range(uc_nat):
        pos_sc_id += [i]*n_uc
    pos_sc_id = np.array(pos_sc_id)
    return pos_sc_id

def dm_calc(q, ph, poscar_sc, poscar_uc, pos_sc_idx): 
    smallest_vectors,multiplicity=ph.primitive.get_smallest_vectors()
    species = poscar_uc.species
    natom = len(species)
    ph.produce_force_constants()
    fc = ph.force_constants
    d_matrix = np.zeros([natom*3, natom*3], dtype = complex)
    for x in range(natom):
        for y in range(natom):
            idx_x = pos_sc_idx[pos_sc_idx==x]
            idx_y = pos_sc_idx[pos_sc_idx==y]
            id1 = 3*x
            id2 = 3*y
            m1 = species[x].atomic_mass
            m2 = species[y].atomic_mass
            for a in range(len(idx_y)): # iterate over all the sublattice y in the supercell
                fc_here = fc[x*len(idx_x),y*len(idx_y)+a]
                multi = multiplicity[y*len(idx_y)+a][x]
                for b in range(multi):
                    vec = smallest_vectors[y*len(idx_y)+a][x][b]
                    vec = vec[0] * A0[:,0] + vec[1] * A0[:,1]
                    exp_here = np.exp(1j * np.dot(q, vec[0:2]))
                    d_matrix[id1:id1+3, id2:id2+3] += fc_here * exp_here / np.sqrt(m1) / np.sqrt(m2) / multi
    d_matrix = (d_matrix + d_matrix.conj().transpose()) / 2 # impose Hermiticity
    return d_matrix
    
main = "/Users/jonathanlu/Documents/tmos2/layer_1/analyses/phonon/"

os.chdir(main)
ph1 = phonopy.load(unitcell_filename="POSCAR_unit", 
                   supercell_matrix=[[3,0,0],[0,3,0],[0,0,1]], 
                   primitive_matrix=np.eye(3),
                   log_level=1)
# help(phonopy.load)
# breakpoint()

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
    id1 = kidx*(nq-1)+1
    for i in range(2):
        kline_tmp = np.linspace(corners[kidx,i], corners[kidx+1,i], nq)
        if kidx != np.shape(corners)[0]-1:
            k_set[id1:id1+nq-1, i] = kline_tmp[0:nq-1]
        else:
            k_set[id1:id1+nq, i] = kline_tmp
    dk = corners[kidx+1,:] - corners[kidx,:]
    dk_norm = norm(dk)
    kmag_tmp = np.linspace(kmag0, kmag0+dk_norm, nq)
    kmag0 = kmag_tmp[-1]
    if kidx != np.shape(corners)[0]-1:
        k_mags[id1:id1+nq-1] = kmag_tmp[0:nq-1]
    else:
        k_mags[id1:id1+nq] = kmag_tmp

corner_kmags = k_mags[0]
for kidx in range(np.shape(corners)[0]-1):
    corner_kmags=np.append([corner_kmags], k_mags[(kidx+1)*(nq-1)])
# fix, ax = plt.subplots()
# ax.scatter(k_set[:,0], k_set[:,1])
# G0_set = np.array([[0,0], G0[:,0], G0[:,1], -G0[:,0], -G0[:,1], G0[:,0]+G0[:,1], -G0[:,1]-G0[:,0]])
# for i, G in enumerate(G0_set[1:]):
#     plt.text(G[0], G[1], str(i+1))
# plt.title("Sampling space"); plt.xlabel(r"$k_x$"); plt.ylabel(r"$k_y$")
# ax.scatter(G0_set[:,0], G0_set[:,1])
# ax.set_aspect("equal")
# plt.show()

# test the functions (work fine)
pos_sc_idx = sc_idx(mos2_sc, mos2)
d_matrix = np.zeros([9, 9, k_set.shape[0]], dtype=complex)
d_matrix2 = np.zeros_like(d_matrix)

evals = np.zeros([9, k_set.shape[0]], dtype=complex)
evals2 = np.zeros_like(evals)

k_direct = np.zeros_like(k_set)

for q_idx, k in enumerate(k_set): 
    k_direct[q_idx] = np.linalg.inv(G0).dot(k)
    
for q_idx, k in enumerate(k_set):
    d_matrix[:, :, q_idx] = dm_calc(k, ph1, mos2_sc, mos2, pos_sc_idx)
    evals[:, q_idx] = 15.633302**2 * np.sort(eigvals(d_matrix[:,:,q_idx]))
    d_matrix2[:, :, q_idx] = ph1.get_dynamical_matrix_at_q(k_direct[q_idx])
    evals2[:, q_idx] = 15.633302**2 * np.real(np.sort(eigvals(d_matrix2[:,:,q_idx])))

fig,ax=plt.subplots(1,2, figsize=(20,6))
ax[0].plot(k_mags, np.transpose(np.sqrt(evals)), color='black')
ax[0].set_title("dm_calc")
ax[0].set_xticks(corner_kmags)
ax[0].set_xticklabels(xticklabels)
ax[0].set_xlim([0, np.max(corner_kmags)])
ax[1].plot(k_mags, np.transpose(np.sqrt(evals)), color='black')
ax[1].set_title("get_dynamical_matrix_at_q")
ax[1].set_xticks(corner_kmags)
ax[1].set_xticklabels(xticklabels)
ax[1].set_xlim([0, np.max(corner_kmags)])
plt.savefig("/Users/jonathanlu/Documents/ml-ph/ALLEGRO_ANALYZER/work_cmp.png")


