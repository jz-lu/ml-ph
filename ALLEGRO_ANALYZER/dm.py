import numpy as np
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
from pymatgen.io.vasp.inputs import Poscar
from itertools import product as prod
from bzsampler import BZSampler

"""
Computes the twisted dynamical matrix from a set of sampled moire G vectors and k points along IBZ boundary.
The large size of the matrix encourages breaking it into blocks. Define a level-0 block to be the intralayer dynamical 
matrix for a specific choice of (GM, k) in a monolayer. A level-1 intralayer block is the concatenation of all such 
level-0 matrices, indexed by the 2 indices of GM (for a given k). A level-0 interlayer block is the sum of over 
all configurations (b-vectors) for a given GM (interlayer terms are independent of k by approximation). A level-1 
interlayer block is exactly analogous to a level-1 intralayer block. Finally, a level-2 block combines the level-1
blocks in the following fashion
      | L1_intra_1     L1_inter     |
L2 =  |                             |    where L{i}_intra_{j} is level i for layer j
      | L1_inter_dag   L1_intra_2   |
There is one level-2 matrix for every k, which are the twisted Fourier dynamical matrices. Diagonalization of each of these
yields the phonon modes as a function of k.
"""
class TDMCalculator:
    def __init__(self, poscar_uc : Poscar, poscar_sc : Poscar, GM_set, k_set):
        self.uc = poscar_uc.structure; self.sc = poscar_sc.structure # get structure objects from Poscar objects
        self.GM_set = GM_set; self.k_set = k_set
        self.pos_sc_idx = self.__sc_idx()
        self.A0 = self.uc.lattice.matrix[:2, :2] # remove z-axis

    # Return supercell index 
    def __sc_idx(self):
        # TODO ??? 90% of the function isn't even used?
        uc_coords = self.uc.cart_coords; sc_coords = self.sc.cart_coords

        sc_nat = len(sc_coords); magic = int(len(sc_coords)/3)
        assert sc_nat % 3 == 0
        pos_m = sc_coords[:magic] - uc_coords[0]; msz = len(pos_m)
        pos_x1 = sc_coords[magic : 2*magic] - uc_coords[1]
        pos_x2 = sc_coords[2*magic:] - uc_coords[2]

        pos_sc_idx = np.zeros(sc_nat)
        for i in range(msz):
            pos_sc_idx[i] = 0
            pos_sc_idx[i + msz] = 1
            pos_sc_idx[i + 2*msz] = 2 # sublattice index
        return pos_sc_idx
    
    # Compute intralayer dynamical matrix block element for given some center `q` and phonopy object `ph`
    def __block_intra_l0(self, q, ph):
        assert self.pos_sc_idx is not None, "Fatal error in class initialization"
        smallest_vectors, multiplicity = ph.primitive.get_smallest_vectors()
        species = self.uc.species; uc_nat = len(species); d = 3 # Cartesian DOF
        fc = ph.force_constants
        D = np.zeros([uc_nat*d, uc_nat*d], dtype=complex) # size: n_at x d (n_at of uc)

        """
        Iteration: choose a pair of atoms (may be the same), fix one, then iterate through all instances
        of the other. TODO ??? 
        """
        for x, y in prod(range(uc_nat), range(uc_nat)):
            idxs_x = self.pos_sc_idx[self.pos_sc_idx == x]; n_xidxs = len(idxs_x)
            idxs_y = self.pos_sc_idx[self.pos_sc_idx == y]; n_yidxs = len(idxs_y)
            id1 = d * x; id2 = d * y
            Mx = species[x].atomic_mass; My = species[y].atomic_mass
            for a in range(n_yidxs): # iterate over all the sublattice y in the supercell
                fc_here = fc[x*n_xidxs, y*n_yidxs + a]
                multi = multiplicity[y*n_yidxs + a][x]
                for vec in smallest_vectors[y*n_yidxs + a][x]:
                    vec = vec[0] * self.A0[:,0] + vec[1] * self.A0[:,1] # convert to Cartesian coords to get R vector
                    fourier_exp = np.exp(1j * np.dot(q, vec[:2]))
                    D[id1:id1+d, id2:id2+d] += (1/np.sqrt(Mx * My)) * fc_here * fourier_exp / multi
        D = (D + D.conj().T) / 2 # impose Hermiticity if not already satisfied
        return D
        
    # Create level-1 block matrix (each matrix becomes a block in level-2)
    def __block_intra_l1(self):
        pass # TODO

    def __block_inter_l0(self):
        pass # TODO

    def __block_inter_l1(self):
        pass # TODO

    # Create level-2 block matrix with intralayer and interlayer terms
    def __block_l2(self, D_intras, D_inter):
        assert len(D_intras) == 2
        assert D_intras[0].shape == D_intras[1].shape == D_inter.shape
        return np.block([[D_intras[0], D_inter], [D_inter.conjudate().T, D_intras[1]]])

    """
    Build the set of dynamical matrices, one for each k in the k-set. Only the intralayer blocks 
    change from one to the next, not the interlayer terms. Returns a list of pairs {k, D(k)}.
    """
    def build_DM_set(self):
        pass
