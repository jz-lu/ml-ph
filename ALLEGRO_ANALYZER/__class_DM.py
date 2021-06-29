import numpy as np
import numpy.linalg as LA
import os
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
from pymatgen.io.vasp.inputs import Poscar
from itertools import product as prod
from bzsampler import BZSampler
import matplotlib.pyplot as plt
from __directory_searchers import checkPath
from ___constants_phonopy import SUPER_DIM
from ___constants_names import DEFAULT_PH_BAND_PLOT_NAME
from scipy.linalg import block_diag
from scipy.sparse import bmat # block matrix

"""
These classes together compute the twisted dynamical matrix from a sample of moire G vectors and k points along IBZ boundary.
The large size of the matrix encourages breaking it into blocks. Define a level-0 block to be the intralayer dynamical 
matrix for a specific choice of (GM, k) in a monolayer. A level-1 intralayer block is the concatenation of all such 
level-0 matrices, indexed by the 2 indices of GM (for a given k). A level-0 interlayer block is the sum of over 
all configurations (b-vectors) for a given GM (interlayer terms are independent of k by approximation). A level-1 
interlayer block is exactly analogous to a level-1 intralayer block. Finally, a level-2 block combines the level-1
blocks in the following fashion

      |   L1_intra_1       L1_inter_12  |
L2 =  |                                 |    where L{i}_intra_{j} is level i for layer j (dag is adjoint)
      | L1_inter_dag_12    L1_intra_2   |

There is one level-2 matrix for every k, which are the twisted Fourier dynamical matrices. Diagonalization of each of these
yields the phonon modes as a function of k.
"""


# Monolayer dynamical matrix
class MonolayerDM:
    def __init__(self, poscar_uc : Poscar, poscar_sc : Poscar, ph, GM_set, k_set, using_flex=False):
        self.uc = poscar_uc.structure; self.sc = poscar_sc.structure # get structure objects from Poscar objects
        self.GM_set = GM_set; self.k_set = k_set; self.ph = ph
        self.pos_sc_id = self.__sc_atomic_id() if using_flex else None
        self.A0 = self.uc.lattice.matrix[:2, :2] # remove z-axis
        self.DM_set = None; self.dbgprint = True
        self.name = poscar_uc.comment
        if self.pos_sc_id is not None:
            print(f"MonolayerDM intralayer atomic IDs for {self.name}:", self.pos_sc_id)

    # Assign each atom a unique ID in the supercell
    def __sc_atomic_id(self):
        uc_coords = self.uc.cart_coords; sc_coords = self.sc.cart_coords
        sc_nat = len(sc_coords); uc_nat = len(uc_coords) # num sc/uc atoms
        pos_sc_id = []; n_uc = int(sc_nat/uc_nat) # num unit cells

        # SPOSCAR arranges all atoms of each type contiguously, so the indexes must
        # be the same for each contiguous region of `n_uc`.
        for i in range(uc_nat):
            pos_sc_id += [i]*n_uc
        pos_sc_id = np.array(pos_sc_id); self.pos_sc_id = pos_sc_id
        return pos_sc_id
    
    # Compute intralayer dynamical matrix block element for given some center `q` and phonopy object `ph`
    def __block_intra_l0(self, q, ph):
        if self.dbgprint:
            print(f"Intralayer Level-0 shape: {ph.get_dynamical_matrix_at_q(q).shape}")
            self.dbgprint = False
        return ph.get_dynamical_matrix_at_q(q)
    
    # Deprecated: Calculates dynamical matrix elements for direct or Cartesian coordinates `q`
    def __flex_block_intra_l0(self, q, ph):
        assert self.pos_sc_id is not None, "Fatal error in class initialization"
        smallest_vectors, multiplicity = ph.primitive.get_smallest_vectors()
        species = self.uc.species; uc_nat = len(species); d = 3 # Cartesian DOF
        fc = ph.force_constants
        D = np.zeros([uc_nat*d, uc_nat*d], dtype=complex) # size: n_at x d (n_at of uc)

        """
        Iteration: choose a pair of atoms (given by sc ID, may be the same), fix one, 
        then iterate through all instances of the other. Add to the matrix the Fourier term
        divided by the multiplicity, which is determined by phonopy based on nearest-neighbor symmetry.
        """
        for x, y in prod(range(uc_nat), range(uc_nat)):
            idxs_x = self.pos_sc_id[self.pos_sc_id == x]; n_xidxs = len(idxs_x)
            idxs_y = self.pos_sc_id[self.pos_sc_id == y]; n_yidxs = len(idxs_y)
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
        self.DM_set = [block_diag(*[self.__block_intra_l0(k+GM, self.ph) for GM in self.GM_set]) for k in self.k_set]
        return self.DM_set

    def get_DM_set(self):
        if self.DM_set is None:
            self.__block_intra_l1()
        print(f"Retrieved monolayer dynamical matrix of shape {self.DM_set[0].shape} for solid {self.name}")
        return self.DM_set

    def get_GM_set(self):
        print(f"Retrieved GM sample set for solid {self.name}")
        return self.GM_set

    def get_k_set(self):
        print(f"Retrieved k sample set for solid {self.name}")
        return self.k_set


# Build interlayer dynamical matrix block via summing over configurations
class InterlayerDM:
    def __init__(self, b_set, ph_list, GM_set, G0_set):
        assert len(b_set[0]) == 2, "Shift vectors must be 2-dimensional"
        self.b_set = b_set; self.ph_list = ph_list # list of phonopy objects for each config
        self.nshift = len(b_set)
        self.GM_set = GM_set; self.G0_set = G0_set
        self.DM = None
        # The force constants matrix for each configuration is equivalent to the dynamical matrix
        # at the Gamma point, since the Fourier transform cancels for G = Gamma.
        self.force_matrices = [ph.get_dynamical_matrix_at_q([0,0,0]) for ph in ph_list]
        assert self.force_matrices[0].shape[0] == self.force_matrices[0].shape[1], f"Force matrix is not square: shape {self.force_matrices[0].shape}"
        assert self.force_matrices[0].shape[0] % 2 == 0, f"Force matrix size is odd: shape {self.force_matrices[0].shape}"
        self.half_pt = self.force_matrices[0].shape[0] // 2

    def __block_inter_l0(self, G0):
        D = sum([force_matrix * np.exp(1j * np.dot(G0, b)) for force_matrix, b in zip(self.force_matrices, self.b_set)])
        return D[:-self.half_pt,-self.half_pt:] / (self.nshift**2) # keep top right corner only

    def __block_inter_l1(self):
        n_GM = len(self.GM_set); assert len(self.G0_set) == n_GM, f"|G0_set| {len(self.G0_set)} != |GM_set| = {n_GM}"
        assert LA.norm(self.G0_set[0]) == 0, f"G0[0] should be 0, but is {LA.norm(self.GM_set[0])}"
        D0 = self.__block_inter_l0(self.G0_set[0]); block_l0_shape = D0.shape
        self.DM = [[None]*n_GM for _ in range(n_GM)] # NoneType interpreted by scipy as 0 matrix block
        for i in range(n_GM): # fill diagonal
            self.DM[i][i] = D0
        for i in range(1, n_GM): # fill first row/col
            self.DM[0][i] = self.__block_inter_l0(self.G0_set[i])
            self.DM[i][0] = self.__block_inter_l0(-self.G0_set[i])
            assert self.DM[0][i].shape == block_l0_shape and self.DM[i][0].shape == block_l0_shape, f"Shape GM0{i}={self.DM[0][i].shape}, GM{i}0={self.DM[i][0].shape}, expected {block_l0_shape}"
            assert np.isclose(LA.norm(self.DM[0][i]), LA.norm(self.DM[i][0]), rtol=1e-5), f"Level-0 interlayer DM blocks for G0{i} not inversion-symmetric:\n {LA.norm(self.DM[0][i])}\nvs. \n{LA.norm(self.DM[i][0])}"
        self.DM = bmat(self.DM).toarray() # convert NoneTypes to zero-matrix blocks to make sparse matrix
        return self.DM

    def get_DM(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        print(f"Retrieved interlayer dynamical matrix of shape {self.DM.shape}")
        return self.DM

    def get_GM_set(self):
        print("Retrieved GM sample set")
        return self.GM_set


# Build full dynamical matrix from intralayer and interlayer terms via the above 2 classes
class TwistedDM:
    def __init__(self, l1 : MonolayerDM, l2 : MonolayerDM, inter : InterlayerDM, k_mags):
        print("Building dynamical matrix intra(er) blocks...")
        DMs_layer1 = l1.get_DM_set(); DMs_layer2 = l2.get_DM_set(); DM_inter = inter.get_DM()
        print("Blocks built.")
        self.DMs = [self.__block_l2([DMs_layer1[i], DMs_layer2[i]], DM_inter) for i in range(len(k_mags))]
        self.k_mags = k_mags
        self.modes_built = False

    # Create level-2 (final level--full matrix) block matrix with intralayer and interlayer terms
    def __block_l2(self, DM_intras, DM_inter):
        assert len(DM_intras) == 2
        assert DM_intras[0].shape == DM_intras[1].shape == DM_inter.shape
        return np.block([[DM_intras[0], DM_inter], [DM_inter.conjugate().T, DM_intras[1]]])
    
    # Retreieve list dynamical matrices corresponding to the list of sampled k-vectors
    def get_DM_set(self):
        print(f"Retrieved DM set of shape {self.DMs[0].shape} from twisted DM object")
        return self.DMs
    
    def get_k_set(self):
        print("Retrieved k set sample from twisted DM object")
        return self.k_mags

    # Diagonalize the set of DMs to get phonon modes
    def build_modes(self):
        self.mode_set = [0]*len(self.k_mags)
        for i, (k_mag, DM) in enumerate(zip(self.k_mags, self.DMs)):
            evals = LA.eigvals(DM); signs = (-1) * (evals < 0) # hack: pull negative sign out of square root to plot imaginary frequencies
            modes_k = signs * np.sqrt(np.abs(evals)) * (15.633302*33.356) # eV/Angs^2 -> THz ~ 15.633302; THz -> cm^-1 ~ 33.356
            self.mode_set[i] = (k_mag, modes_k)
        # self.k_plt = []; self.modes_plt = []
        # for k_mag, modes in self.mode_set:
        #     self.k_plt += [k_mag]*len(modes)
        #     self.modes_plt += list(modes)
        self.modes_built = True
        return self.mode_set
    
    # Plot phonon modes as a function of k
    def plot_band(self, corner_kmags, angle, outdir='./', filename=DEFAULT_PH_BAND_PLOT_NAME, name=None):
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        if not self.modes_built:
            print("Modes not built yet, building...")
            self.build_modes()
            print("Modes built.")
        plt.clf()
        for k_mag, modes in self.mode_set:
            plt.scatter([k_mag] * len(modes), modes, c='blue')
        xlabs = (r'$\Gamma$', r'K', r'M', r'$\Gamma$')
        plt.xlabel(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = r"Phonon modes"
        if name is not None:
            title += f" of {name} bilayer"
        title += r" at " + '%.1lf'%angle + r"$^\circ$"
        plt.title(title)
        plt.savefig(outdir + filename)
        print(f"Plotting completed and written to {outdir+filename}")
        return


