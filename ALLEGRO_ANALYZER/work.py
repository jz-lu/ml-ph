import phonopy, os
from pymatgen.io.vasp.inputs import Poscar
import numpy as np
import numpy.linalg as LA
from math import pi
import matplotlib.pyplot as plt

def sc_idx(poscar_sc, poscar_uc):
    pos_sc = poscar_sc.cart_coords
    pos_uc = poscar_uc.cart_coords
    A0 = np.transpose(poscar_uc.lattice.matrix[0:2,0:2])

    pos_m = pos_sc[0:int(len(pos_sc)/3)]-pos_uc[0,:]
    pos_x1 = pos_sc[int(len(pos_sc)/3):2*int(len(pos_sc)/3)]-pos_uc[1,:]
    pos_x2 = pos_sc[2*int(len(pos_sc)/3):]-pos_uc[2,:]
    pos_sc_idx = np.zeros([len(pos_sc)])
    
    for i in range(len(pos_m)):
        pos_sc_idx[i] = 0
        pos_sc_idx[i+len(pos_m)] = 1
        pos_sc_idx[i+2*len(pos_m)] = 2 # sublattice index
    return pos_sc_idx

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

a = '/Users/jonathanlu/Documents/tmos2/layer_1/analyses/phonon/'
uc = Poscar.from_file(a+'POSCAR_unit')
sc = Poscar.from_file(a+'SPOSCAR')
A0 = uc.structure.lattice.matrix[:2,:2].T
G0 = 2 * pi * LA.inv(A0).T
os.chdir(a)
ph = phonopy.load(supercell_filename='SPOSCAR')
FREQUENCY_CONVERSTION_FACTOR = 521.4644215120001
d = 2

nq = nk = 101
K = [0.666667, 0.333333]
M = [0.5, 0.5]
Gamma = [0, 0]
pt = corners = np.array([K, Gamma, M, K]); ncorners = len(pt)

nsample = (nk-1) * (ncorners-1)
kline = np.zeros([nsample, d]); kmag = np.zeros(nsample); kmag_start = 0
corner_kmags = []
for line in range(ncorners-1): # last point equals first, so skip it
    kidx = line*(nk-1) # convert line index to k-index
    # Drop second corner point in each line to avoid resampling corner points
    kline[kidx : kidx+nk-1] = np.linspace(corners[line], corners[line+1], nk)[:-1]
    dline_mag = LA.norm(corners[line+1] - corners[line])
    mags = np.linspace(kmag_start, kmag_start + dline_mag, nk)
    corner_kmags.append(kmag_start)
    kmag_start = mags[-1] # update start point of flattened-k to end of current line
    kmag[kidx : kidx+nk-1] = mags[:-1]
Gamma_idx = nk-1

pos_sc_idx = sc_idx(sc.structure, uc.structure)
# evals = np.array([np.real_if_close(LA.eigvals(dm_calc(q, ph, sc.structure, uc.structure, pos_sc_idx))) for q in kline])
evals = np.array([np.real_if_close(LA.eigvals(ph.get_dynamical_matrix_at_q(q))) for q in kline])
evals = (-1*(evals < 0) + (evals > 0)) * np.sqrt(np.abs(evals))
evals = 15.633302 * evals
for q, modes in zip(kmag, evals):
    plt.scatter([q]*len(modes), modes, s=0.1, color='black')
plt.xticks(corner_kmags, ["K", r"$\Gamma$", "M"])
plt.ylabel("Frequency (THz)")
plt.savefig('/Users/jonathanlu/Documents/ml-ph/ALLEGRO_ANALYZER/work.png')

