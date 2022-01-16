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
from ___constants_sampling import DEFAULT_KDIM, DEFAULT_PARTITION_DENSITY, DOS_NEGLIGIBLITY_PROP
from ___constants_phonopy import SUPER_DIM
from ___constants_names import DEFAULT_PH_BAND_PLOT_NAME, DEFAULT_PH_BANDDOS_PLOT_NAME
from ___constants_vasp import VASP_FREQ_TO_INVCM_UNITS
from ___helpers_parsing import update, succ
from scipy.linalg import block_diag
from scipy.sparse import bmat # block matrix
from scipy.signal import savgol_filter as smoothen_filter
from math import sqrt, pi, log, ceil
from copy import deepcopy as dc
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

def sc_idx(sc, uc):
    uc_coords = uc.cart_coords; sc_coords = sc.cart_coords
    sc_nat = len(sc_coords); uc_nat = len(uc_coords) # num sc/uc atoms
    pos_sc_id = []; n_uc = int(sc_nat/uc_nat) # num unit cells

    # SPOSCAR arranges all atoms of each type contiguously, so the indices must
    # be the same for each contiguous region of `n_uc`.
    for i in range(uc_nat):
        pos_sc_id += [i]*n_uc
    pos_sc_id = np.array(pos_sc_id)
    return pos_sc_id

def dm_calc(q, ph, poscar_sc, poscar_uc, pos_sc_idx): 
    smallest_vectors,multiplicity=ph.primitive.get_smallest_vectors()
    species = poscar_uc.species
    A0 = poscar_uc.lattice.matrix[:2,:2].T
    natom = len(species)
    fc = ph.force_constants
    if fc.shape[0] != fc.shape[1]:
        ph.produce_force_constants()
        fc = ph.force_constants
        assert fc.shape[0] == fc.shape[1]
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
        # q = LA.inv(self.G0) @ q
        # dm = ph.get_dynamical_matrix_at_q(q)
        dm = dm_calc(q, ph, self.sc, self.uc, sc_idx(self.sc, self.uc))
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
            blks = [dm_calc([0,0], self.ph, self.sc, self.uc, sc_idx(self.sc, self.uc)) for _ in self.G0_set]
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
        plt.clf(); bot = 100
        for k_mag, modes in self.mode_set:
            if cutoff is not None:
                modes = modes[modes <= cutoff]
            bot = min(bot, min(modes))
            plt.scatter([k_mag] * len(modes), modes, c='royalblue', s=0.07)
        xlabs = (r'K', r'$\Gamma$', r'M', r'K')
        plt.xticks(corner_kmags, xlabs)
        plt.ylim(bottom=bot)
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
                 species_per_layer, sc, uc, 
                 ph_list=None, 
                 force_matrices=None, G0=None, br_set=None, A0=None):
        self.DM_G_blks = None
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
        self.bl_M = bl_M; self.bl_n_at = len(bl_M)
        assert self.bl_n_at == len(uc.species)
        self.ph_list = ph_list
        self.br_set = br_set # relaxed b vector set

        print(f"Interlayer DM uses bilayer with {self.bl_n_at} atoms of masses {self.bl_M}")

        assert (ph_list is not None) ^ (force_matrices is not None), "Must give exactly one of: phonopy obj list, force matrix list"
        assert len(b_set[0]) == 2, "Shift vectors must be 2-dimensional"
        self.b_set = b_set; self.ph_list = ph_list # list of phonopy objects for each config
        self.nshift = len(b_set); print(f"Number of configurations: {self.nshift}")
        assert int(sqrt(self.nshift))**2 == self.nshift, f"Number of shifts {self.nshift} must be a perfect square"
        self.GM_set = GM_set; self.G0_set = G0_set; self.DM = None; self.DM_intra = None
        self.per_layer_at_idxs = per_layer_at_idxs; assert len(self.per_layer_at_idxs) == 2, f"Only 2 layers supported"
        cat_pl_idxs = np.concatenate(self.per_layer_at_idxs)
        self.k_set = k_set
        self.modes_built = False
        self.n_GM = len(GM_set)
        if ph_list is not None:
            pos_sc_idx = sc_idx(sc, uc)
            self.force_matrices = [dm_calc([0,0], ph, sc, uc, pos_sc_idx) for ph in ph_list]
            # self.force_matrices = [ph.get_dynamical_matrix_at_q([0,0]) for ph in ph_list] # DM(Gamma) ~= FC (mass-scaled)
        else:
            self.force_matrices = force_matrices
        self.force_matrices = [fm[cat_pl_idxs] for fm in self.force_matrices]
        self.force_matrices = [fm[:,cat_pl_idxs] for fm in self.force_matrices]
        def symmetrize_intralayer(fm):
            mid = fm.shape[0] // 2
            fm[:mid, :mid] = fm[mid:, mid:] = 1/2 * (fm[:mid, :mid] + fm[mid:, mid:])
            return fm
        # self.force_matrices = [symmetrize_intralayer(fm) for fm in self.force_matrices]
        
        if br_set is not None:
            assert A0 is not None, f"Must give lattice matrix A0 when using relaxer"
            br_sz = int(sqrt(self.br_set.shape[0]))
            update("CONFIGURATION DM: Using RELAXED configs")
            relaxed_forces = [self.__inverse_block_inter_l0(br) for br in self.br_set]
            self.force_matrices = relaxed_forces
            a = np.linspace(0,1,br_sz, endpoint=False)
            self.b_set = np.array(list(prod(a,a))) @ A0.T
       
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
        cut = len(self.per_layer_at_idxs[0])
        D = sum([force_matrix * np.exp(1j * np.dot(G0 + k, b)) \
                        for force_matrix, b in zip(self.force_matrices, self.b_set)])
        # Extract a submatrix with rows of atoms from layer 1 
        # and columns of atoms from layer 2, which is the interlayer 1-2 interactions.
        # D_inter = D[np.ix_(self.per_layer_at_idxs[0], self.per_layer_at_idxs[1])] / self.nshift
        # D_intra1 = D[np.ix_(self.per_layer_at_idxs[0], self.per_layer_at_idxs[0])] / self.nshift
        # D_intra2 = D[np.ix_(self.per_layer_at_idxs[1], self.per_layer_at_idxs[1])] / self.nshift
        D_inter = D[:cut,cut:] / self.nshift
        D_intra1 = D[:cut,:cut] / self.nshift
        D_intra2 = D[cut:,cut:] / self.nshift
        self.l0sz = D_inter.shape
        assert len(D.shape) == 2 and D.shape[0] == D.shape[1], f"D with shape {D.shape} not square matrix"
        assert len(D_inter.shape) == 2 and D_inter.shape[0] == D_inter.shape[1], f"D_inter with shape {D_inter.shape} not square matrix"
        assert 2*D_inter.shape[0] == D.shape[0] and 2*D_inter.shape[1] == D.shape[1], f"D_inter shape {D_inter.shape} should be half of D shape {D.shape}"
        return D_inter, D_intra1, D_intra2
    
    def __inverse_block_inter_l0(self, br):
        if self.DM_G_blks is None:
            # Get the blocks over all G0
            self.DM_G_blks = [sum([force_matrix * np.exp(1j * np.dot(G0, b)) \
                        for force_matrix, b in zip(self.force_matrices, self.b_set)]) / self.nshift \
                        for G0 in self.G0_set]
        # Fourier transform into the given br
        br = br[:2]; assert br.shape == (2,), f"Invalid relaxed b shape {br.shape}"
        assert len(self.DM_G_blks) == len(self.G0_set), \
            f"len(self.DM_G_blks) = {len(self.DM_G_blks)} != len(self.G0_set) = {len(self.G0_set)}"
        return sum([dm_G * np.exp(-1j * np.dot(G0, br)) \
                        for dm_G, G0 in zip(self.DM_G_blks, self.G0_set)])

    def __block_inter_l1(self):
        n_GM = self.n_GM; assert len(self.G0_set) == n_GM, f"|G0_set| {len(self.G0_set)} != |GM_set| = {n_GM}"
        assert LA.norm(self.G0_set[0]) == 0, f"G0[0] should be 0, but is {LA.norm(self.GM_set[0])}"
        D0, D0intra1, D0intra2 = self.__block_inter_l0(self.G0_set[0]); block_l0_shape = D0.shape
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
            # self.DM_intra[0][i][i] = D0intra1
            # self.DM_intra[1][i][i] = D0intra2
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
        
        # DM_inter = np.zeros_like(DM_inter); DM_cfgintra = np.zeros_like(DM_cfgintra) #! delete!!
        
        self.szs = [DMs_layer1[0].shape[0], DMs_layer2[0].shape[0]]
        print("Blocks built.")
        
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
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
        self.sum_rule_applied = False
        self.intra_set = [block_diag(DMs_layer1[i], DMs_layer2[i]) for i in range(self.n_k)]
        assert self.intra_set[0].shape == self.DMs[0].shape
    
    # Dynamical matrix without any interlayer coupling
    def get_intra_set(self): 
        return np.array(self.intra_set)

    # Create level-2 (final level--full matrix) block matrix with intralayer and interlayer terms
    def __block_l2(self, DM_intras, DM_inter):
        assert len(DM_intras) == 2
        assert DM_intras[0].shape == DM_intras[1].shape == DM_inter.shape
        DM = np.block([[DM_intras[0], DM_inter], [DM_inter.conjugate().T, DM_intras[1]]])
        return DM

    def apply_sum_rule(self):
        update("ENFORCING: acoustic sum rule")
        if self.sum_rule_applied:
            update("Already applied sum rule!")
            return
        else:
            self.sum_rule_applied = True
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

    def get_DM_set(self):
        """Retrieve list dynamical matrices corresponding to the list of sampled k-vectors"""
        print(f"Retrieved DM set of shape {self.DMs[0].shape} from twisted DM object")
        return np.array(self.DMs)
    
    def get_DM_at_Gamma(self):
        """Dynamical matrix at Gamma point"""
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

    def build_modes(self, dump=False, outdir='.'):
        """Diagonalize the set of DMs to get phonon modes"""
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
    
    def get_mode_set(self):
        if not self.modes_built:
            self.build_modes()
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
                assert len(modes) > 0, f"Cutoff size {cutoff} too small!"
            lims[0] = min(0, lims[0], min(modes))
            lims[1] = max(lims[1], max(modes))
            plt.scatter([k_mag] * len(modes), modes, c='black', s=0.07)
        lims[1] += 5; plt.ylim(lims)
        xlabs = (r'K', r'$\Gamma$', r'M', r'K')
        plt.xticks(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = r"Phonon modes"
        if name is not None:
            title += f" of {name}"
        title += r" at " + '%.3lf'%angle + r"$^\circ$"
        plt.title(title)
        plt.savefig(outdir + filename, format='pdf')
        plt.close()
        print(f"Band plot written to {outdir+filename}")
        return

"""Analysis of dynamical matrix for a given k-point in theta-space"""
class ThetaSpaceDM:
    def __init__(self, k, TDMs, thetas, modeidxs=np.linspace(0,6,6)):
        assert len(thetas) == len(TDMs), \
            f"Number of angles {len(thetas)} does not math number of DMs {len(TDMs)}"
        self.DM_set = TDMs; self.thetas = thetas
        self.modeidxs = modeidxs
        print(f"Initialized dynamical matrix analyzer in theta-space for k = {k}.")
        self.__analyze()
        
    def __DM_to_DDM(self):
        """Diagonalize and sort eigensystem"""
        def sorted_filtered_eigsys(A):
            evals, evecs = LA.eig(A)
            evals = np.real(evals); evecs = np.real_if_close(evecs)
            idxs = evals.argsort()  
            evals = evals[idxs]; evecs = evecs[:,idxs].T
            evals = evals[self.modeidxs]; evecs = evecs[self.modeidxs]
            return evals, evecs
        
        # Let C = len(modeidxs) and n_th = len(thetas)
        eigsys = [sorted_filtered_eigsys(A) for A in self.DM_set]
        evals = np.array([vals for vals,_ in eigsys]) # shape: (n_th, C)
        self.DDM = np.transpose([vecs for _,vecs in eigsys], axes=(1,2,0)) # shape: (C, n_at x n_G x d, n_th)
        self.modes = np.sign(evals) * np.sqrt(np.abs(evals)) * VASP_FREQ_TO_INVCM_UNITS
        print("Diagonalization complete.")
    
    def __DDM_to_thtnsr(self):
        """Diagonalized dynamical matrix -> tensor object. Run after __DM_to_DDM"""
        evecs_by_layer = np.split(self.DDM, 2, axis=1) # layer slice must be first
        for li in range(2): # each layer has shape: (C, n_G x n_at/2 x d, n_th)
            evecs_by_layer[li] = np.array(np.split(evecs_by_layer[li], self.n_G, axis=1)) # shape: (n_G, C, n_at/2 x d, n_th)
            evecs_by_layer[li] = np.array(np.split(evecs_by_layer[li], self.n_at//2, axis=2)) # shape: (n_at/2, n_G, C, d, n_th)
            print(f"Final evec shape for layer {li}: {evecs_by_layer[li].shape}")
        evecs = np.concatenate(evecs_by_layer, axis=0) # bring them together, so shape: (n_at, n_G, C, d, n_th)
        self.thtnsr = np.transpose(evecs, axes=(1,0,2,3,4)) # shape: (n_G, n_at, C, d, n_th)
        print("Tensors prepared.")

    def __analyze(self):
        self.__DM_to_DDM()
        self.__DDM_to_thtnsr()
        print("Analysis in theta-space complete.")
    
    def get_modes(self):
        return self.modes

    def get_thtnsr(self):
        return self.thtnsr
        

class TwistedDOS:
    def __init__(self, TDMs, n_G, theta, cutoff=None, 
                 partition_density=DEFAULT_PARTITION_DENSITY, 
                 kdim=DEFAULT_KDIM, normalizer=1, eigsys=None):
        # Eigensort the eigenbasis, then slice off everything but the G0 component
        def sorted_sliced_eigsys(A, idx, n_G):
            vals, vecs = LA.eig(A); idxs = vals.argsort()   
            vals = vals[idxs]; vecs = vecs[:,idxs]
            mid = vecs.shape[0] // 2; glen = vecs.shape[0] // (n_G * 2)
            # Cut matrix in half (split by layer), then take the G0 (First) component
            vecs = np.vstack((vecs[0 : 0+glen], vecs[mid : mid+glen]))
            assert vecs.shape == (glen*2, vals.shape[0])
            if idx % 100 == 0:
                print(".", end='', flush=True)
            return vals, vecs

        assert TDMs[0].shape[0] % (n_G * 2) == 0
        self.theta = theta
        print(f"Initialized Twisted DOS for theta={self.theta}")

        self.eigsys = eigsys
        if self.eigsys is None:
            print("Diagonalizing Moire dynamical matrices", end='', flush=True)
            self.eigsys = [sorted_sliced_eigsys(TDM, i, n_G) for i, TDM in enumerate(TDMs)]
            print(" done", flush=True)
        self.modes = np.array([(-1*(evals < 0) + 1*(evals > 0)) \
            * np.sqrt(np.abs(evals)) * VASP_FREQ_TO_INVCM_UNITS for evals,_ in self.eigsys])
        self.weights = np.array([np.square(LA.norm(evecs, axis=0)) for _,evecs in self.eigsys])
        assert self.modes.shape == self.weights.shape == (kdim**2, TDMs[0].shape[0])
        self.modes = self.modes.flatten(); self.weights = self.weights.flatten()
        self.mode_extrema = [np.min(self.modes), np.max(self.modes)]
        if cutoff is not None:
            print(f"DOS CUTOFF: {cutoff}")
            assert cutoff > self.mode_extrema[0], \
                f"Cutoff {cutoff} smaller than minimum freq {self.mode_extrema[0]}"
            self.mode_extrema[1] = min(self.mode_extrema[1], cutoff)
        self.__set_parameters(theta, partition_density, normalizer)
        self.DOS = None
    
    def get_eigsys(self):
        return self.eigsys
    
    def __obtain_width(self, theta):
        if theta < 1.01:
            return 0.03
        if theta < 1.51:
            return 0.09
        if theta < 2.01:
            return 0.13
        if theta < 2.51:
            return 0.15
        if theta < 3.51:
            return 0.21
        if theta < 5.51:
            return 0.27
        if theta < 6.51:
            return 0.31
        if theta < 7.51:
            return 0.37
        if theta < 9.01:
            return 0.43
        else:
            return 0.45
        
    def __set_parameters(self, theta, partition_density, normalizer):
        # Gaussian is f(x) = A e^(x^2/2 sigma)
        # Width parameter converts to sigma via http://hyperphysics.phy-astr.gsu.edu/hbase/Math/gaufcn2.html 
        self.bin_sz = (self.mode_extrema[1] - self.mode_extrema[0]) / partition_density
        self.omegas = np.linspace(self.mode_extrema[0], self.mode_extrema[1], num=partition_density)
        self.npts = partition_density
        self.width = self.__obtain_width(theta)
        self.sigma = self.width / (2 * sqrt(2 * log(2))) 
        self.A = normalizer
        print(f"DOS PARAMETERS:\n\tA: {self.A}\n\tdE: {self.bin_sz}\n\tWidth: {self.width}\n\tsigma: {self.sigma}")

    # Compute DOS at some list of omegas
    def __DOS_at_omegas(self, omegas):
        def DOS_at_omega(enum, A, sigma, modes, weights, bin_sz):
            idx, omega = enum
            if idx % 100 == 0:
                print(".", end="", flush=True)
            return sum([weight * A * np.exp(-np.power(omega - omega_k, 2.) / (2 * np.power(sigma, 2.))) \
                for omega_k, weight in zip(modes, weights) ])

        print("Getting DOS at omegas", end="", flush=True)
        DOSs = np.array([DOS_at_omega(enum, self.A, self.sigma, \
            self.modes, self.weights, self.bin_sz) for enum in enumerate(omegas)])
        print(" done", flush=True)
        return DOSs

    # Obtain DOS for plotting / any other use
    def get_DOS(self, smoothen=True, wsz=33, polyd=3):
        if self.DOS is None:
            self.DOS = self.__DOS_at_omegas(self.omegas)
            self.max_DOS = np.max(self.DOS)
        if smoothen:
            smooth_DOS = dc(self.DOS)
            scan_wsz = self.npts // 40
            nwindow = int(ceil(len(self.DOS) / scan_wsz))
            for i in range(nwindow):
                start, end = scan_wsz*i, scan_wsz*(i+1)
                window = self.DOS[start:end]
                wmin, wmax = np.min(window), np.max(window)
                if abs(wmax-wmin) <= self.max_DOS * DOS_NEGLIGIBLITY_PROP:
                    wsz_here = min(wsz, len(window))
                    if wsz_here % 2 == 0:
                        wsz_here -= 1
                    if wsz_here <= 0:
                        continue
                    smooth_DOS[start:end] = smoothen_filter(window, wsz_here, polyd)
            # smooth_DOS = smoothen_filter(self.DOS, wsz, polyd) # window size 31, polynomial order 7
            return self.omegas, smooth_DOS
        else:
            return self.omegas, self.DOS
    
    def plot_DOS(self, vertical=True, outdir='.'):
        assert os.path.isdir(outdir), f"Invalid directory '{outdir}'"
        outpath = checkPath(os.path.abspath(outdir)) + f"dos{self.theta}_{self.width}.pdf"
        plt.clf(); fig, ax = plt.subplots()
        self.get_DOS()
        if vertical:
            ax.plot(self.DOS, self.omegas, c='black', linewidth=0.9)
            ax.set_xlabel("DOS")
            ax.set_ylabel(r"$\omega (cm^{-1})$")
        else:
            ax.plot(self.omegas, self.DOS, c='black', linewidth=0.9)
            ax.set_ylabel("DOS")
            ax.set_xlabel(r"$\omega$ $(cm^{-1})$")
        ax.set_title(rf"DOS ($\theta = {self.theta}$" \
            + r"$^\circ, width = $" + f"{self.width})")
        fig.savefig(outpath)
        succ(f"Successfully outputted DOS-{self.theta} to {outpath}")
        return self.omegas, self.DOS, ax
    

class TwistedPlotter:
    def __init__(self, theta, omegas, DOS, mode_set, corner_kmags, cutoff=None):
        self.theta = theta
        self.omegas = omegas
        self.DOS = DOS
        self.mode_set = mode_set
        self.corner_kmags = corner_kmags
        self.cutoff = cutoff
        print("Initialized TwistedPlotter object")
    
    def make_plot(self, outdir='./', filename=DEFAULT_PH_BANDDOS_PLOT_NAME, name=None, width=None, pfx=None):
        if width is not None:
            filename = "w%f_"%(round(width, 3)) + filename
        if pfx is not None:
            filename = pfx + "_" + filename
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        plt.clf(); fig, [axband, axdos] = plt.subplots(nrows=1, ncols=2, sharey=True)
        lims = [0,0]
        for k_mag, modes in self.mode_set:
            if self.cutoff is not None:
                modes = modes[modes <= self.cutoff]
            lims[0] = min(0, lims[0], min(modes))
            lims[1] = max(lims[1], max(modes))
            axband.scatter([k_mag] * len(modes), modes, c='black', s=0.07)
        lims[1] += 0.05*lims[1]; plt.ylim(lims)
        xlabs = (r'K', r'$\Gamma$', r'M', r'K')
        axband.set_xticks(self.corner_kmags)
        axband.set_xticklabels(xlabs)
        axband.set_ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = name + " b" if name is not None else "B"
        title += "and " + '(%.1lf'%self.theta + r"$^\circ$)"
        axband.set_title(title)

        if self.cutoff is not None:
            self.DOS = self.DOS[self.omegas <= self.cutoff]
            self.omegas = self.omegas[self.omegas <= self.cutoff]
        axdos.plot(self.DOS, self.omegas, c='black', linewidth=0.6)
        # axdos.scatter(self.DOS, self.omegas, c='black', s=1)
        axdos.set_title("DOS")
        axdos.ticklabel_format(axis='x', scilimits=(0,0), style='sci')
        plt.ylim(lims)
        plt.savefig(outdir + filename, format='pdf')
        plt.close(fig)
        print(f"Band-DOS plot written to {outdir+filename}")

        