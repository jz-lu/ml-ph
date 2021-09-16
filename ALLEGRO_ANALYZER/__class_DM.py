import numpy as np
import numpy.linalg as LA
import os
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
from ___helpers_parsing import update, succ
from scipy.linalg import block_diag
from scipy.sparse import bmat # block matrix
from math import sqrt, pi
import sys
from __class_PhonopyAPI import PhonopyAPI

"""
These classes tosgether compute the twisted dynamical matrix from a sample of moire G vectors and k points along IBZ boundary.
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
    def __init__(self, poscar_uc : Poscar, poscar_sc : Poscar, ph, 
                 GM_set, G0_set, k_set, Gamma_idx, k_mags=None):
        assert LA.norm(k_set[Gamma_idx]) == 0, f"Gamma index has nonzero norm {LA.norm(k_set[Gamma_idx])}"
        print("Symmetrizing force constants...")
        ph.symmetrize_force_constants()
        print("Force constants symmetrized.")
        self.n_at = sum(poscar_uc.natoms)
        self.uc = poscar_uc.structure; self.sc = None if poscar_sc is None else poscar_sc.structure # get structure objects from Poscar objects
        self.GM_set = GM_set; self.n_GM = len(GM_set); self.k_set = np.array(k_set); self.ph = ph
        self.G0_set = G0_set
        assert np.all(k_set[Gamma_idx] == [0,0]), f"Nonzero kGamma={k_set[Gamma_idx]}"
        self.n_k = len(k_set)
        self.A0 = self.uc.lattice.matrix[:2,:2].T # remove z-axis
        self.G0 = 2 * pi * LA.inv(self.A0).T
        self.DM_set = None; self.dbgprint = True; self.l0_shape = None
        self.name = poscar_uc.comment; self.modes_built = False
        self.M = np.array([species.atomic_mass for species in poscar_uc.structure.species])
        self.Gamma_idx = Gamma_idx
        self.k_mags = k_mags
        self.l0_shape = [3*self.n_at]*2
        self.corr_mat = None
        # super().__init__(self.sc, self.uc)
    
    # Compute intralayer dynamical matrix block element for given some (direct) center `q` and phonopy object `ph`
    def __block_intra_l0(self, q, ph):
        q = LA.inv(self.G0) @ q
        dm = ph.get_dynamical_matrix_at_q(q)
        if self.dbgprint:
            print(f"Intralayer Level-0 shape: {dm.shape}")
            self.dbgprint = False
        return dm
    
    # Create level-1 block matrix (each matrix becomes a block in level-2)
    def __block_intra_l1(self):
        self.DM_set = [block_diag(*[self.__block_intra_l0(k+GM, self.ph) for GM in self.GM_set]) for k in self.k_set]
        return self.DM_set
    
    def Gamma_blks(self, log=False):
        blks = [self.__block_intra_l0(GM, self.ph) for GM in self.GM_set]
        if log:
            for i, blk in enumerate(blks):
                print(f"G{i}", np.real_if_close(blk))
        return blks

    def get_block_force_sum(self):
        if self.DM_set is None:
            self.__block_intra_l1()
        l0sz = self.l0_shape[0]; Gam = self.Gamma_idx
        def force_from_dm(i):
            Mi = self.M[i]
            temp = np.split(self.DM_set[Gam][3*i:3*(i+1),:l0sz], self.n_at, axis=1)
            return sum([sqrt(Mi*Mj) * subblk for subblk, Mj in zip(temp, self.M)])
        blksums_list = [np.real_if_close(force_from_dm(i)) for i in range(self.n_at)]
        return blksums_list
    
    def print_force_sum(self):
        blksums_list = self.get_block_force_sum()
        print(f"Masses: {self.M}")
        print(f"\nINTRALAYER FORCE SUMS for {self.name}:")
        blksum_str = [f"[At: {i}]:\n{s}\n" for i, s in enumerate(blksums_list)]
        print(''.join(blksum_str))

    def get_DM_set(self):
        if self.DM_set is None:
            self.__block_intra_l1()
        print(f"Retrieved monolayer dynamical matrix of shape {self.DM_set[0].shape} for solid {self.name}")
        print(f"Frobenius norm of pristine intralayer piece at Gamma: {LA.norm(self.DM_set[self.Gamma_idx])}")
        return self.DM_set
    
    def get_corr_mat(self):
        if self.corr_mat is None:
            blks = [self.ph.get_dynamical_matrix_at_q([0,0]) for _ in self.G0_set]
            self.corr_mat = block_diag(*blks)
        return self.corr_mat
    
    def get_origin_G_blocks(self):
        if self.DM_set is None:
            self.__block_intra_l1()
        # blks = [self.__block_intra_l0(GM, self.ph) for GM in self.GM_set if LA.norm(GM) > 0]
        # assert np.allclose(block_diag(self.__block_intra_l0(0, self.ph), *blks), self.DM_set[self.Gamma_idx])
        l0sz = self.l0_shape[0]
        DM = self.DM_set[self.Gamma_idx]
        blks = [DM[i*l0sz:(i+1)*l0sz, i*l0sz:(i+1)*l0sz] for i in range(1, self.n_at)]
        print(f"Retrieved {len(blks)} blocks from intralayer G-blocks")
        return blks
    
    def plot_sampled_l0(self, outdir='.'):
        outdir = checkPath(os.path.abspath(outdir))
        lst = np.absolute([self.__block_intra_l0(GM, self.ph)[0,0] for GM in self.GM_set])
        plt.clf(); fig, ax = plt.subplots()
        ax.scatter(np.arange(self.n_GM), lst)
        fig.savefig(outdir + 'sampledl0.png')
        plt.close(fig)
    
    def plot_pristine_band(self, k0_set, k0_mags, corner_k0mags, outdir='.'):
        d_matrix = np.zeros([9, 9, k0_set.shape[0]], dtype=complex)
        evals = np.zeros([9, k0_set.shape[0]], dtype=complex)
        
        for q_idx, k0 in enumerate(k0_set):
            d_matrix[:, :, q_idx] = self.__block_intra_l0(k0 + self.G0_set[1], self.ph)
            evals[:, q_idx] = VASP_FREQ_TO_INVCM_UNITS**2 * np.real(np.sort(LA.eigvals(d_matrix[:,:,q_idx])))

        _,ax = plt.subplots(1,1, figsize=(8,6))
        evals = (-1*(evals < 0) + 1*(evals > 0)) * np.sqrt(np.abs(evals))
        ax.plot(k0_mags, np.transpose(evals), color='black')
        ax.set_title(f"ph from API")
        ax.set_xticks(corner_k0mags)
        ax.set_xticklabels(["K", r"$\Gamma$", "M", "K"])
        ax.set_xlim([0, np.max(corner_k0mags)])
        plt.savefig(outdir + "pristine.png")
    
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
            self.mode_set[i] = (k_mag, modes_k)
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
        xlabs = (r'K', r'$\Gamma$', r'M', r'K')
        plt.xticks(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = "Moire intralayer phonon modes"
        if name is not None:
            title += f" of {name}"
        plt.title(title)
        plt.savefig(outdir + filename)
        succ(f"Moire intralayer plot written to {outdir+filename}")
        return

    def get_GM_set(self):
        print(f"Retrieved GM sample set for solid {self.name}")
        return self.GM_set

    def get_k_set(self):
        print(f"Retrieved k sample set for solid {self.name}")
        return self.k_set


# Build interlayer dynamical matrix block via summing over configurations
class InterlayerDM:
    def __init__(self, per_layer_at_idxs, bl_M, 
                 b_set, k_set, 
                 GM_set, G0_set, 
                 species_per_layer, 
                 ph_list=None, 
                 force_matrices=None, G0=None):
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
        self.bl_M = bl_M; self.bl_n_at = len(bl_M)
        self.ph_list = ph_list
        print(f"Interlayer DM uses bilayer with {self.bl_n_at} atoms of masses {self.bl_M}")

        assert (ph_list is not None) ^ (force_matrices is not None), "Must give exactly one of: phonopy obj list, force matrix list"
        assert len(b_set[0]) == 2, "Shift vectors must be 2-dimensional"
        self.b_set = b_set; self.ph_list = ph_list # list of phonopy objects for each config
        self.nshift = len(b_set); print(f"Number of configurations: {self.nshift}")
        assert int(sqrt(self.nshift))**2 == self.nshift, f"Number of shifts {self.nshift} must be a perfect square"
        self.GM_set = GM_set; self.G0_set = G0_set; self.DM = None; self.DM_intra = None
        self.per_layer_at_idxs = per_layer_at_idxs; assert len(self.per_layer_at_idxs) == 2, f"Only 2 layers supported"
        self.k_set = k_set
        self.modes_built = False
        self.n_GM = len(GM_set)
        if ph_list is not None:
            self.force_matrices = [ph.get_dynamical_matrix_at_q([0,0]) for ph in ph_list] # DM(Gamma) ~= FC (mass-scaled)
        else:
            self.force_matrices = force_matrices
       
        assert self.force_matrices[0].shape[0] == self.force_matrices[0].shape[1], f"Force matrix is not square: shape {self.force_matrices[0].shape}"
        assert self.force_matrices[0].shape[0] % 2 == 0, f"Force matrix size is odd: shape {self.force_matrices[0].shape}"

    # Plot pristine bilayer for a given shift (currently shift=30 which is AB stacking)
    def plot_pristine_band(self, G0_mat, k0_set, k0_mags, corner_k0mags, outdir='.'):
        modes = np.zeros((18,len(k0_set)))
        for kidx, k in enumerate(k0_set):
            dm = self.ph_list[0].get_dynamical_matrix_at_q(LA.inv(G0_mat) @ k)
            evals = np.sort(LA.eigvals(dm))
            modes[:,kidx] = (-1*(evals < 0) + 1*(evals > 0))*np.sqrt(np.abs(evals))
        _,ax = plt.subplots(1,1, figsize=(8,6))
        ax.plot(k0_mags, np.transpose(modes), color='black')
        ax.set_title(f"AA ph inter from API")
        ax.set_xticks(corner_k0mags)
        ax.set_xticklabels(["K", r"$\Gamma$", "M", "K"])
        ax.set_xlim([0, np.max(corner_k0mags)])
        plt.savefig(outdir + "blpristine.png")
        succ("ILDM DONE")

    def __block_inter_l0(self, G0, k=np.array([0,0])):
        D = sum([force_matrix * np.exp(-1j * np.dot(G0 + k, b)) for force_matrix, b in zip(self.force_matrices, self.b_set)])
        # Extract a submatrix with rows of atoms from layer 1 
        # and columns of atoms from layer 2, which is the interlayer 1-2 interactions.
        D_inter = D[np.ix_(self.per_layer_at_idxs[0], self.per_layer_at_idxs[1])] / self.nshift
        D_intra1 = D[np.ix_(self.per_layer_at_idxs[0], self.per_layer_at_idxs[0])] / self.nshift
        D_intra2 = D[np.ix_(self.per_layer_at_idxs[1], self.per_layer_at_idxs[1])] / self.nshift
        self.l0sz = D_inter.shape
        assert len(D.shape) == 2 and D.shape[0] == D.shape[1], f"D with shape {D.shape} not square matrix"
        assert len(D_inter.shape) == 2 and D_inter.shape[0] == D_inter.shape[1], f"D_inter with shape {D_inter.shape} not square matrix"
        assert 2*D_inter.shape[0] == D.shape[0] and 2*D_inter.shape[1] == D.shape[1], f"D_inter shape {D_inter.shape} should be half of D shape {D.shape}"
        return D_inter, D_intra1, D_intra2

    def __block_inter_l1(self):
        n_GM = self.n_GM; assert len(self.G0_set) == n_GM, f"|G0_set| {len(self.G0_set)} != |GM_set| = {n_GM}"
        assert LA.norm(self.G0_set[0]) == 0, f"G0[0] should be 0, but is {LA.norm(self.GM_set[0])}"
        D0,_,_ = self.__block_inter_l0(self.G0_set[0]); block_l0_shape = D0.shape
        # D0 = 1.15 * D0
        self.DM = [[None]*n_GM for _ in range(n_GM)] # NoneType interpreted by scipy as 0 matrix block
        self.DM_intra = [[[None]*n_GM for _ in range(n_GM)] for _ in range(2)]
        self.GMi_intra_blocks = [[0]*(n_GM-1) for i in range(2)]
        self.all_blocks = [0]*(n_GM)
        self.all_blocks[0] = D0
        self.Gamma_block = D0
        for i in range(1, n_GM): # fill first row/col
            self.DM[0][i], self.DM_intra[0][0][i], self.DM_intra[1][0][i] = self.__block_inter_l0(self.G0_set[i])
            self.GMi_intra_blocks[0][i-1] = self.DM_intra[0][0][i]
            self.GMi_intra_blocks[1][i-1] = self.DM_intra[1][0][i]
            self.all_blocks[i] = self.DM[0][i]
            self.DM[i][0], self.DM_intra[0][i][0], self.DM_intra[1][i][0] = self.__block_inter_l0(-self.G0_set[i])
            assert self.DM[0][i].shape == block_l0_shape and self.DM[i][0].shape == block_l0_shape, f"Shape GM0{i}={self.DM[0][i].shape}, GM{i}0={self.DM[i][0].shape}, expected {block_l0_shape}"
            assert np.isclose(LA.norm(self.DM[0][i]), LA.norm(self.DM[i][0]), rtol=1e-5), f"Level-0 interlayer DM blocks for G0{i} not inversion-symmetric:\n {LA.norm(self.DM[0][i])}\nvs. \n{LA.norm(self.DM[i][0])}"
            for j in range(2):
                assert np.isclose(LA.norm(self.DM_intra[j][0][i]), LA.norm(self.DM_intra[j][i][0]), rtol=1e-5), \
                    f"Level-0 cfg-intra-l{j} DM blocks for G0{i} not inversion-symmetric:\
                        \n {LA.norm(self.DM_intra[j][0][i])}\nvs. \n{LA.norm(self.DM_intra[j][i][0])}"
        for i in range(n_GM): # fill diagonal, but only for interlayer piece
            self.DM[i][i] = D0
        self.DM = bmat(self.DM).toarray() # convert NoneTypes to zero-matrix blocks to make sparse matrix
        for i in range(2):
            self.DM_intra[i] = bmat(self.DM_intra[i]).toarray()
            print(f"Frobenius norm of cfg-{i} intralayer piece at Gamma: {LA.norm(self.DM_intra[i])}")
        return self.DM

    # Diagonalize the set of DMs to get phonon modes
    def __build_modes(self, k_mags, DM_set, dump=False, outdir=None):
        self.k_mags = k_mags
        self.DM_set = DM_set
        self.mode_set = [0]*len(self.k_mags)
        for i, (k_mag, DM) in enumerate(zip(self.k_mags, self.DM_set)):
            evals = LA.eigvals(DM)
            signs = (-1)*(evals < 0) + (evals > 0) # pull negative sign out of square root to plot imaginary frequencies
            modes_k = signs * np.sqrt(np.abs(evals)) * (VASP_FREQ_TO_INVCM_UNITS) # eV/Angs^2 -> THz ~ 15.633302; THz -> cm^-1 ~ 33.356
            self.mode_set[i] = (k_mag, modes_k)
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
    
    def get_off_diag_blocks(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        print(f"Off-diagonal level-0 size: {self.l0sz}")
        self.off_diag_blks = np.split(self.DM[:self.l0sz[0]], self.n_GM, axis=1)[1:]
        print(f"Off-diagonal block ({len(self.off_diag_blks)} ct.) shapes: {self.off_diag_blks[0].shape}")
        return self.off_diag_blks
    
    def get_Gamma_block(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        print(f"Interlayer Gamma block shape: {self.Gamma_block.shape}")
        return self.Gamma_block
    
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
    
    def get_DM_intra_piece(self):
        if self.DM is None:
            print("Building interlayer dynamical matrix...")
            self.__block_inter_l1()
        print(f"Retrieved interlayer dynamical matrix of shape {self.DM.shape}")
        assert self.DM_intra[0].shape == self.DM_intra[1].shape == self.DM.shape
        return self.DM_intra

    def get_GM_set(self):
        print("Retrieved GM sample set")
        return self.GM_set


# Build full dynamical matrix from intralayer and interlayer terms via the above 2 classes
class TwistedDM:
    def __init__(self, l1 : MonolayerDM, l2 : MonolayerDM, inter : InterlayerDM, k_mags, species_per_layer, Gamma_idx):
        self.interobj = inter; self.intraobjs = [l1, l2]
        self.n_ats = [l1.n_at, l2.n_at]
        self.GM_set = inter.get_GM_set()
        print("Building dynamical matrix intra(er) blocks...")
        DMs_layer1 = l1.get_DM_set(); DMs_layer2 = l2.get_DM_set()
        DM_inter = inter.get_DM(); DM_cfgintra = inter.get_DM_intra_piece()
        # DMs_layer1, DMs_layer2 = inter.get_intra_DM_set()
        self.szs = [DMs_layer1[0].shape[0], DMs_layer2[0].shape[0]]
        print("Blocks built.")
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
        # self.off_diag_blocks = inter.get_all_blocks()
        self.off_diag_blocks = inter.get_off_diag_blocks()
        self.on_diag_blocks = [l.get_origin_G_blocks() for l in [l1, l2]]
        self.intra_config_blocks = inter.get_intra_blocks()
        self.n_GM = len(l1.get_GM_set())
        self.k_mags = k_mags; self.n_k = len(k_mags)
        self.modes_built = False
        self.l0szs = [l1.l0_shape[0], l2.l0_shape[0]]
        self.Gamma_idx = Gamma_idx
        self.Gamma_intra_blks = [l1.Gamma_blks(), l2.Gamma_blks()]
        self.DMs = [self.__block_l2([DMs_layer1[i] + DM_cfgintra[0],\
             DMs_layer2[i] + DM_cfgintra[1]], DM_inter) for i in range(self.n_k)]
        self.corr_mat = self.__block_l2([l1.get_corr_mat() + DM_cfgintra[0], \
             l2.get_corr_mat() + DM_cfgintra[1]], DM_inter)

    # Create level-2 (final level--full matrix) block matrix with intralayer and interlayer terms
    def __block_l2(self, DM_intras, DM_inter):
        assert len(DM_intras) == 2
        assert DM_intras[0].shape == DM_intras[1].shape == DM_inter.shape
        DM = np.block([[DM_intras[0], DM_inter], [DM_inter.conjugate().T, DM_intras[1]]])
        return DM

    def apply_sum_rule(self):
        update("ENFORCING: Lukas sum rule")
        dm = self.corr_mat # similar to D(Gamma) but intralayer terms are all at Gamma instead of GM
        blkshp = len(dm)//3
        corr = np.zeros_like(dm, dtype=complex)
        M = self.M.tolist()
        M = M[0]*self.n_GM + M[1]*self.n_GM # masses per block
        assert len(M) == blkshp, f"Inconsistent Cartesian-reduced moire DM shape and mass list shapes {blkshp} vs. {len(M)}"
        for i in range(blkshp):
            corr[3*i:3*(i+1),3*i:3*(i+1)] = corr[3*i:3*(i+1),3*i:3*(i+1)] - sum([sqrt(M[j]/M[i]) \
                * dm[3*j:3*(j+1), 3*i:3*(i+1)] \
                    for j in range(blkshp)])
        for k in range(self.n_k):
            self.DMs[k] = self.DMs[k] + corr

    # Retrieve list dynamical matrices corresponding to the list of sampled k-vectors
    def get_DM_set(self):
        print(f"Retrieved DM set of shape {self.DMs[0].shape} from twisted DM object")
        return np.array(self.DMs)
    
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
        interblks = [self.interobj.get_Gamma_block()]
        print(f"Using {len(interblks)} blocks from interlayer component")
        print("LAYER 1 ATOMS:")
        for i in range(self.n_ats[0]): # layer 1 atoms
            intra = dm[:self.szs[0],:self.szs[0]]
            intrasum = np.real_if_close(intra_force_from_dm(0, i, 0, l0szs[0], self.n_ats[0], intra))
            assert intrasum.shape == (3,3)
            intersum = np.real_if_close(sum([inter_force_from_dm(0, i, self.n_ats[1], blk) for blk in interblks]))
            assert intersum.shape == (3,3)
            totalsum = np.real_if_close(np.round(intrasum + intersum, 8))
            print(f"[Layer 1] [At {i}] total sum of reciprocal-forces at Gamma:\n{totalsum}")
            # print(f"[Layer 1] [At {i}] intra sum of forces:\n{totalsum}")
            
        print("Layer 2 ATOMS:")
        for i in range(self.n_ats[1]): # layer 2 atoms
            intra = dm[self.szs[0]:self.szs[0]+self.szs[1],self.szs[0]:self.szs[0]+self.szs[1]]
            assert dm.shape == tuple(map(lambda x: 2*x, intra.shape))
            intrasum = intra_force_from_dm(1, i, 0, l0szs[1], self.n_ats[1], intra)
            intersum = sum([inter_force_from_dm(1, i, self.n_ats[0], blk) for blk in interblks])
            assert intrasum.shape == (3,3), f"Invalid intra shape {intrasum.shape}"
            assert intersum.shape == (3,3), f"Invalid inter shape {intersum.shape}"
            totalsum = np.real_if_close(np.round(intrasum + intersum, 8))
            print(f"[Layer 2] [At {i}] total sum of reciprocal-forces at Gamma:\n{totalsum}")

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
            self.mode_set[i] = (k_mag, modes_k)
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
        lims = [0,0]
        for k_mag, modes in self.mode_set:
            if cutoff is not None:
                modes = modes[modes <= cutoff]
            lims[0] = min(0, lims[0], min(modes))
            lims[1] = max(lims[1], max(modes))
            plt.scatter([k_mag] * len(modes), modes, c='black', s=0.07)
        lims[1] += 30; plt.ylim(lims)
        xlabs = (r'K', r'$\Gamma$', r'M', r'K')
        plt.xticks(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = r"Phonon modes"
        if name is not None:
            title += f" of {name} bilayer"
        title += r" at " + '%.1lf'%angle + r"$^\circ$"
        plt.title(title)
        plt.savefig(outdir + filename, format='pdf')
        print(f"Plot written to {outdir+filename}")
        return


