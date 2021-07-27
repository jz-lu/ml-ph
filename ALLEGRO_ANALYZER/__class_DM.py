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
from ___constants_vasp import VASP_FREQ_TO_INVCM_UNITS
from scipy.linalg import block_diag
from scipy.sparse import bmat # block matrix
from math import sqrt

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
    def __init__(self, poscar_uc : Poscar, poscar_sc : Poscar, ph, GM_set, k_set, Gamma_idx, using_flex=False):
        assert LA.norm(k_set[Gamma_idx]) == 0, f"Gamma index has nonzero norm {LA.norm(k_set[Gamma_idx])}"
        print("Symmetrizing force constants...")
        ph.symmetrize_force_constants()
        print("Force constants symmetrized.")
        self.n_at = sum(poscar_uc.natoms)
        self.uc = poscar_uc.structure; self.sc = poscar_sc.structure # get structure objects from Poscar objects
        self.GM_set = GM_set; self.n_GM = len(GM_set); self.k_set = k_set; self.ph = ph
        self.pos_sc_id = self.__sc_atomic_id() if using_flex else None
        self.A0 = self.uc.lattice.matrix[:2, :2] # remove z-axis
        self.DM_set = None; self.dbgprint = True; self.l0_shape = None
        self.name = poscar_uc.comment; self.modes_built = False
        self.M = np.array([species.atomic_mass for species in poscar_uc.structure.species])
        if self.pos_sc_id is not None:
            print(f"MonolayerDM intralayer atomic IDs for {self.name}:", self.pos_sc_id)
        self.Gamma_idx = Gamma_idx

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
    
    # Compute intralayer dynamical matrix block element for given some (direct) center `q` and phonopy object `ph`
    def __block_intra_l0(self, q, ph):
        dm = ph.get_dynamical_matrix_at_q(q)
        def enforce_trnsl_invrnc():
            if LA.norm(q) == 0.0:
                print(f"Translational invariance ignored for q={q}")
                return
            assert dm.shape[0] == dm.shape[1] and dm.shape[0] == 3*self.n_at, f"Shape {dm.shape} inconsistent with num atoms {self.n_at}"
            for i in range(0, dm.shape[0], 3):
                dm[i:i+3,i:i+3] = np.zeros_like(dm[i:i+3,i:i+3])
                    
        if self.dbgprint:
            print(f"Intralayer Level-0 shape: {dm.shape}")
            self.dbgprint = False
        if self.l0_shape is None:
            self.l0_shape = dm.shape
        enforce_trnsl_invrnc()
        return dm
    
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
    
    def Gamma_blks(self):
        return [self.__block_intra_l0(GM, self.ph) for GM in self.GM_set]

    def get_block_force_sum(self):
        if self.DM_set is None:
            self.__block_intra_l1()
        l0sz = self.l0_shape[0]; Gam = self.Gamma_idx
        def force_from_dm(i, GMj):
            Mi = self.M[i]; idx = 3*i + GMj*l0sz
            temp = np.split(self.DM_set[Gam][idx:idx+3,GMj*l0sz:(GMj+1)*l0sz], self.n_at, axis=1)
            return sum([sqrt(Mi*Mj) * subblk for subblk, Mj in zip(temp, self.M)])
        if self.DM_set is None:
            self.__block_intra_l1()
        blksums_list = [np.round(force_from_dm(i, 0), 8) for i in range(self.n_at)]
        return blksums_list
    
    def print_force_sum(self):
        blksums_list = self.get_block_force_sum()
        print(f"\nINTRALAYER FORCE SUMS for {self.name}:")
        blksum_str = [f"[At: {i}]:\n{s}\n" for i, s in enumerate(blksums_list)]
        print(''.join(blksum_str))

    def get_DM_set(self):
        if self.DM_set is None:
            self.__block_intra_l1()
        print(f"Retrieved monolayer dynamical matrix of shape {self.DM_set[0].shape} for solid {self.name}")
        return self.DM_set
    
    # Diagonalize the set of DMs to get phonon modes
    def build_modes(self, k_mags, dump=False, outdir=None):
        self.k_mags = k_mags
        if self.DM_set is None:
            self.__block_intra_l1()
        self.mode_set = [0]*len(self.k_mags)
        for i, (k_mag, DM) in enumerate(zip(self.k_mags, self.DM_set)):
            evals = LA.eigvals(DM)
            signs = (-1)*(evals < 0) + (evals > 0) # pull negative sign out of square root to plot imaginary frequencies
            modes_k = signs * np.sqrt(np.abs(evals)) * (VASP_FREQ_TO_INVCM_UNITS) # eV/Angs^2 -> THz ~ 15.633302; THz -> cm^-1 ~ 33.356
            self.mode_set[i] = (k_mag, modes_k[modes_k != 0])
        if dump:
            outdir = checkPath(os.path.abspath(outdir))
            assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
            k_dump = []; mode_dump = []
            for k_mag, modes in self.mode_set:
                    k_dump += [k_mag]*len(modes)
                    mode_dump += list(modes)
            with open(outdir + 'intra_modes.txt', 'w') as f:
                for k_mag, mode in zip(k_dump, mode_dump):
                    f.write(f"{k_mag}\t{mode}\n")
        self.modes_built = True
        return self.mode_set
    
    def plot_band(self, k_mags, corner_kmags, outdir='./', 
                  filename='intra_' + DEFAULT_PH_BAND_PLOT_NAME, name=None, cutoff=None):
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        if not self.modes_built:
            print("Intralayer modes not built yet, building...")
            self.build_modes(k_mags, dump=False, outdir=outdir)
            print("Intralayer modes built.")
        plt.clf()
        for k_mag, modes in self.mode_set:
            if cutoff is not None:
                modes = modes[modes <= cutoff]
            plt.scatter([k_mag] * len(modes), modes, c='royalblue', s=0.07)
        xlabs = (r'$\Gamma$', r'K', r'M')
        plt.xticks(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = "First-layer phonon modes"
        if name is not None:
            title += f" of {name}"
        plt.title(title)
        plt.savefig(outdir + filename)
        print(f"Intralayer pplot written to {outdir+filename}")
        return

    def get_GM_set(self):
        print(f"Retrieved GM sample set for solid {self.name}")
        return self.GM_set

    def get_k_set(self):
        print(f"Retrieved k sample set for solid {self.name}")
        return self.k_set


# Build interlayer dynamical matrix block via summing over configurations
class InterlayerDM:
    def __init__(self, per_layer_at_idxs, b_set, GM_set, G0_set, species_per_layer, ph_list=None, force_matrices=None):
        assert (ph_list is not None) ^ (force_matrices is not None), "Must give exactly one of: phonopy obj list, force matrix list"
        assert len(b_set[0]) == 2, "Shift vectors must be 2-dimensional"
        self.b_set = b_set; self.ph_list = ph_list # list of phonopy objects for each config
        self.nshift = len(b_set); print(f"Number of configurations: {self.nshift}")
        assert int(sqrt(self.nshift))**2 == self.nshift, f"Number of shifts {self.nshift} must be a perfect square"
        self.GM_set = GM_set; self.G0_set = G0_set; self.DM = None
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
        self.per_layer_at_idxs = per_layer_at_idxs; assert len(self.per_layer_at_idxs) == 2, f"Only 2 layers supported"
        if ph_list is not None:
            def sym_dm_at_gamma(i, ph):
                # print(f"Symmetrizing force constants for config {i}...")
                ph.symmetrize_force_constants()
                return ph.get_dynamical_matrix_at_q([0,0,0])
            self.force_matrices = [sym_dm_at_gamma(i, ph) for i, ph in enumerate(ph_list)] # DM(Gamma) = FC (mass-scaled)
        else:
            self.force_matrices = force_matrices
        assert self.force_matrices[0].shape[0] == self.force_matrices[0].shape[1], f"Force matrix is not square: shape {self.force_matrices[0].shape}"
        assert self.force_matrices[0].shape[0] % 2 == 0, f"Force matrix size is odd: shape {self.force_matrices[0].shape}"

    def __block_inter_l0(self, G0, k=np.array([0,0,0])):
        D = sum([force_matrix * np.exp(1j * np.dot(G0 + k, b)) for force_matrix, b in zip(self.force_matrices, self.b_set)])
        # Extract a submatrix with rows of atoms from layer 1 
        # and columns of atoms from layer 2, which is the interlayer 1-2 interactions.
        D_inter = D[np.ix_(self.per_layer_at_idxs[0], self.per_layer_at_idxs[1])] / self.nshift
        D_intra1 = D[np.ix_(self.per_layer_at_idxs[0], self.per_layer_at_idxs[0])] / self.nshift
        D_intra2 = D[np.ix_(self.per_layer_at_idxs[1], self.per_layer_at_idxs[1])] / self.nshift
        assert len(D.shape) == 2 and D.shape[0] == D.shape[1], f"D with shape {D.shape} not square matrix"
        assert len(D_inter.shape) == 2 and D_inter.shape[0] == D_inter.shape[1], f"D_inter with shape {D_inter.shape} not square matrix"
        assert 2*D_inter.shape[0] == D.shape[0] and 2*D_inter.shape[1] == D.shape[1], f"D_inter shape {D_inter.shape} should be half of D shape {D.shape}"
        return D_inter, D_intra1, D_intra2

    def __block_inter_l1(self):
        # def enforce_acoustic_sum_rule(D):
        #     n_at_1 = len(self.per_layer_at_idxs[0])
        #     assert n_at_1 % 3 == 0
        #     M1 = self.M[0]; M2 = self.M[1]
        #     for i in range(0, 3*n_at_1, 3):
        #         D[i:i+3,i:i+3] -= sum([sum([sqrt(M2[j//3] / M1[i//3]) * Gblk[i:i+3,j:j+3] for j in range(0,Gblk.shape[1],3)]) for Gblk in self.GMi_blocks])
        #     return D

        n_GM = len(self.GM_set); assert len(self.G0_set) == n_GM, f"|G0_set| {len(self.G0_set)} != |GM_set| = {n_GM}"
        assert LA.norm(self.G0_set[0]) == 0, f"G0[0] should be 0, but is {LA.norm(self.GM_set[0])}"
        D0,_,_ = self.__block_inter_l0(self.G0_set[0]); block_l0_shape = D0.shape
        self.DM = [[None]*n_GM for _ in range(n_GM)] # NoneType interpreted by scipy as 0 matrix block
        self.GMi_blocks = [0]*(n_GM-1)
        self.GMi_intra_blocks = [[0]*(n_GM-1) for i in range(2)]
        self.all_blocks = [None]*(3*(n_GM-1) + 1)
        for i in range(1, n_GM): # fill first row/col
            self.DM[0][i], self.GMi_intra_blocks[0][i-1], self.GMi_intra_blocks[1][i-1] = self.__block_inter_l0(self.G0_set[i])
            self.GMi_blocks[i-1] = self.DM[0][i]
            self.all_blocks[i+n_GM-1] = self.DM[0][i]
            self.DM[i][0],_,_ = self.__block_inter_l0(-self.G0_set[i])
            self.all_blocks[i+(2*n_GM-1)-1] = self.DM[i][0]
            assert self.DM[0][i].shape == block_l0_shape and self.DM[i][0].shape == block_l0_shape, f"Shape GM0{i}={self.DM[0][i].shape}, GM{i}0={self.DM[i][0].shape}, expected {block_l0_shape}"
            assert np.isclose(LA.norm(self.DM[0][i]), LA.norm(self.DM[i][0]), rtol=1e-5), f"Level-0 interlayer DM blocks for G0{i} not inversion-symmetric:\n {LA.norm(self.DM[0][i])}\nvs. \n{LA.norm(self.DM[i][0])}"
        
        # D0 = enforce_acoustic_sum_rule(D0)
        for i in range(n_GM): # fill diagonal
            self.DM[i][i] = D0
            self.all_blocks[i] = D0
        self.DM = bmat(self.DM).toarray() # convert NoneTypes to zero-matrix blocks to make sparse matrix
        assert np.round(sum(sum(self.DM)), 12) == np.round(sum(sum(sum(self.all_blocks))), 12), f"The list of blocks does not match the interlayer DM, {sum(sum(sum(self.all_blocks)))} vs. {sum(sum(self.DM))}"
        return self.DM
    
    def get_off_diag_blocks(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        print(f"Off-diagonal block shapes: {self.GMi_blocks[0].shape}")
        return self.GMi_blocks
    
    def get_intra_blocks(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        print(f"Config intralayer block shapes: {self.GMi_intra_blocks[0][0].shape} and {self.GMi_intra_blocks[1][0].shape}")
        return self.GMi_intra_blocks
    
    def get_all_blocks(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        return self.all_blocks

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
    def __init__(self, l1 : MonolayerDM, l2 : MonolayerDM, inter : InterlayerDM, k_mags, species_per_layer, Gamma_idx):
        self.interobj = inter; self.intraobjs = [l1, l2]
        self.n_ats = [l1.n_at, l2.n_at]
        print("Building dynamical matrix intra(er) blocks...")
        DMs_layer1 = l1.get_DM_set(); DMs_layer2 = l2.get_DM_set(); DM_inter = inter.get_DM()
        self.szs = [DMs_layer1[0].shape[0], DMs_layer2[0].shape[0]]
        print("Blocks built.")
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
        self.off_diag_blocks = inter.get_off_diag_blocks()
        self.intra_config_blocks = inter.get_intra_blocks()
        self.n_GM = len(l1.get_GM_set())
        self.k_mags = k_mags
        self.modes_built = False
        self.l0szs = [l1.l0_shape[0], l2.l0_shape[0]]
        self.Gamma_idx = Gamma_idx
        self.Gamma_intra_blks = [l1.Gamma_blks(), l2.Gamma_blks()]
        self.DMs = [self.__block_l2([DMs_layer1[i], DMs_layer2[i]], DM_inter) for i in range(len(k_mags))]

    # Create level-2 (final level--full matrix) block matrix with intralayer and interlayer terms
    def __block_l2(self, DM_intras, DM_inter):
        # * On the first intralayer loop over the 3x3 block diagonal i and sum over all 3x3 blocks
        # * in the interlayer 1-2 matrix in row i. Repeat for second intralayer and interlayer 2-1.
        def enforce_acoustic_sum_rule():
            M = self.M; n_GM = self.n_GM
            l0sz = DM_intras[0].shape[0] // n_GM; assert DM_intras[0].shape[0] % n_GM == 0; assert l0sz % 3 == 0
            dg1 = self.Gamma_intra_blks[0][1:]; dg2 = self.Gamma_intra_blks[1][1:]
            d1 = self.intra_config_blocks[0]; d2 = self.intra_config_blocks[1]
            for i in range(0, l0sz, 3): # intralayer1 and interlayer12
                DM_intras[0][i:i+3,i:i+3] -= sum([sum([sqrt(M[1,j//3] / M[0,i//3]) * Gblk[i:i+3,j:j+3] for j in range(0,Gblk.shape[1],3)]) for Gblk in self.off_diag_blocks])
                # DM_intras[0][i:i+3,i:i+3] -= sum([sum([sqrt(M[0,j//3] / M[0,i//3]) * Gblk[i:i+3,j:j+3] for j in range(0,Gblk.shape[1],3) if j != i]) for Gblk in d1])
                # DM_intras[0][i:i+3,i:i+3] -= sum([sum([sqrt(M[0,j//3] / M[0,i//3]) * Gblk[i:i+3,j:j+3] for j in range(0,Gblk.shape[1],3) if j != i]) for Gblk in dg1])
            l0sz = DM_intras[1].shape[0] // n_GM; assert DM_intras[1].shape[0] % n_GM == 0; assert l0sz % 3 == 0
            for i in range(0, l0sz, 3): # intralayer2 and interlayer 21
                DM_intras[1][i:i+3,i:i+3] -= sum([sum([sqrt(M[0,j//3] / M[1,i//3]) * Gblk.conjugate().T[i:i+3,j:j+3] for j in range(0,Gblk.shape[0],3)]) for Gblk in self.off_diag_blocks])
                # DM_intras[0][i:i+3,i:i+3] -= sum([sum([sqrt(M[1,j//3] / M[1,i//3]) * Gblk[i:i+3,j:j+3] for j in range(0,Gblk.shape[1],3) if j != i]) for Gblk in d2])
                # DM_intras[1][i:i+3,i:i+3] -= sum([sum([sqrt(M[1,j//3] / M[1,i//3]) * Gblk[i:i+3,j:j+3] for j in range(0,Gblk.shape[0],3) if j != i]) for Gblk in dg2])
                
        assert len(DM_intras) == 2
        assert DM_intras[0].shape == DM_intras[1].shape == DM_inter.shape
        DM_intras = list(map(lambda D: (D + D.conjugate().T) / 2, DM_intras)) # impose Hermiticity
        enforce_acoustic_sum_rule()
        DM = np.block([[DM_intras[0], DM_inter], [DM_inter.conjugate().T, DM_intras[1]]])
        return DM
    
    # Retreieve list dynamical matrices corresponding to the list of sampled k-vectors
    def get_DM_set(self):
        print(f"Retrieved DM set of shape {self.DMs[0].shape} from twisted DM object")
        return self.DMs
    
    def get_DM_at_Gamma(self):
        return self.DMs[self.Gamma_idx]
    
    def print_force_sum(self, G0_only=True):
        def intra_force_from_dm(l, i, j, l0sz, n_at, DMintra):
            assert l in [0,1]; Mi = self.M[l,i]
            idx = 3*i + j*l0sz
            temp = np.split(DMintra[idx:idx+3,j*l0sz:(j+1)*l0sz], n_at, axis=1)
            return sum([sqrt(Mi*Mj) * subblk for subblk, Mj in zip(temp, self.M[l])])
        def inter_force_from_dm(l, i, n_at, blk): # n_at is number in *other* layer
            assert l in [0,1]; Mi = self.M[l,i]
            if l == 1:
                blk = blk.conjugate().T
            temp = np.split(blk[3*i:3*(i+1)], n_at, axis=1)
            return sum([sqrt(Mi*Mj) * subblk for subblk, Mj in zip(temp, self.M[1-l])])

        print("\nTWISTED FORCE SUMS:")
        dm = self.get_DM_at_Gamma()
        l0szs = self.l0szs
        interblks = self.interobj.get_all_blocks()
        print(f"Using {len(interblks)} blocks from interlayer component")
        print("LAYER 1 ATOMS:")
        for i in range(self.n_ats[0]): # layer 1 atoms
            intra = dm[:self.szs[0],:self.szs[0]]
            intrasum = intra_force_from_dm(0, i, 0, l0szs[0], self.n_ats[0], intra)
            assert intrasum.shape == (3,3)
            intersum = sum([inter_force_from_dm(0, i, self.n_ats[1], blk) for blk in interblks])
            assert intersum.shape == (3,3)
            totalsum = np.round(intrasum + intersum, 8)
            print(f"[Layer 1] [At {i}] total sum of forces:\n{totalsum}")
            
        print("Layer 2 ATOMS:")
        for i in range(self.n_ats[1]): # layer 2 atoms
            intra = dm[self.szs[0]:self.szs[0]+self.szs[1],self.szs[0]:self.szs[0]+self.szs[1]]
            assert dm.shape == tuple(map(lambda x: 2*x, intra.shape))
            intrasum = intra_force_from_dm(1, i, 0, l0szs[1], self.n_ats[1], intra)
            intersum = sum([inter_force_from_dm(1, i, self.n_ats[0], blk) for blk in interblks])
            assert intrasum.shape == (3,3), f"Invalid intra shape {intrasum.shape}"
            assert intersum.shape == (3,3), f"Invalid inter shape {intersum.shape}"
            totalsum = np.round(intrasum + intersum, 8)
            print(f"[Layer 2] [At {i}] total sum of forces:\n{totalsum}")

    def get_k_set(self):
        print("Retrieved k set sample from twisted DM object")
        return self.k_mags

    # Diagonalize the set of DMs to get phonon modes
    def build_modes(self, dump=False, outdir='.'):
        self.mode_set = [0]*len(self.k_mags)
        self.modetnsr = [0]*len(self.k_mags)
        for i, (k_mag, DM) in enumerate(zip(self.k_mags, self.DMs)):
            evals = LA.eigvals(DM)
            signs = (-1)*(evals < 0) + (evals > 0) # pull negative sign out of square root to plot imaginary frequencies
            modes_k = signs * np.sqrt(np.abs(evals)) * (VASP_FREQ_TO_INVCM_UNITS) # eV/Angs^2 -> THz ~ 15.633302; THz -> cm^-1 ~ 33.356
            self.mode_set[i] = (k_mag, modes_k[modes_k != 0])
            self.modetnsr[i] = modes_k
        if dump:
            outdir = checkPath(os.path.abspath(outdir))
            assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
            k_dump = []; mode_dump = []
            for k_mag, modes in self.mode_set:
                    k_dump += [k_mag]*len(modes)
                    mode_dump += list(modes)
            with open(outdir + 'modes.txt', 'w') as f:
                for k_mag, mode in zip(k_dump, mode_dump):
                    f.write(f"{k_mag}\t{mode}\n")
        self.modes_built = True
        self.modetnsr = np.array(self.modetnsr)
        return self.mode_set
    
    def k_mode_tensor(self):
        if not self.modes_built:
            self.build_modes()
        return self.modetnsr
    
    # Plot phonon modes as a function of k
    def plot_band(self, corner_kmags, angle, outdir='./', filename=DEFAULT_PH_BAND_PLOT_NAME, name=None, cutoff=None):
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        if not self.modes_built:
            print("Modes not built yet, building...")
            self.build_modes(dump=True, outdir=outdir)
            print("Modes built.")
        plt.clf()
        minmax = [0,0]
        for k_mag, modes in self.mode_set:
            if cutoff is not None:
                modes = modes[modes <= cutoff]
            minmax[0] = min(0, minmax[0], min(modes))
            minmax[1] = max(minmax[1], max(modes))
            plt.scatter([k_mag] * len(modes), modes, c='royalblue', s=0.07)
        minmax[1] += 30; plt.ylim(minmax)
        xlabs = (r'K', r'$\Gamma$', r'M')
        plt.xticks(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = r"Phonon modes"
        if name is not None:
            title += f" of {name} bilayer"
        title += r" at " + '%.1lf'%angle + r"$^\circ$"
        plt.title(title)
        plt.savefig(outdir + filename)
        print(f"Plot written to {outdir+filename}")
        return


