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
from ___helpers_parsing import update, succ
from scipy.linalg import block_diag
from scipy.sparse import bmat # block matrix
from math import sqrt, pi
import sys

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

class DM:
    def __init__(self, struct_sc, struct_uc):
        self.sc = struct_sc; self.uc = struct_uc
        self.pos_sc_idx = self.sc_idx()
        print("Initialized parent DM calculator")
    
    def sc_idx(self):
        pos_sc = self.sc.cart_coords
        pos_uc = self.uc.cart_coords
        self.A0 = np.transpose(self.uc.lattice.matrix[0:2,0:2])

        pos_m = pos_sc[0:int(len(pos_sc)/3)]-pos_uc[0,:]
        pos_x1 = pos_sc[int(len(pos_sc)/3):2*int(len(pos_sc)/3)]-pos_uc[1,:]
        pos_x2 = pos_sc[2*int(len(pos_sc)/3):]-pos_uc[2,:]
        pos_sc_idx = np.zeros([len(pos_sc)])
        
        for i in range(len(pos_m)):
            pos_sc_idx[i] = 0
            pos_sc_idx[i+len(pos_m)] = 1
            pos_sc_idx[i+2*len(pos_m)] = 2 # sublattice index
        return pos_sc_idx

    def dm_calc(self, q, ph): 
        pos_sc_idx = self.pos_sc_idx
        smallest_vectors,multiplicity=ph.primitive.get_smallest_vectors()
        species = self.uc.species
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
                        vec = vec[0] * self.A0[:,0] + vec[1] * self.A0[:,1]
                        exp_here = np.exp(1j * np.dot(q, vec[0:2]))
                        d_matrix[id1:id1+3, id2:id2+3] += fc_here * exp_here / np.sqrt(m1) / np.sqrt(m2) / multi
        d_matrix = (d_matrix + d_matrix.conj().transpose()) / 2 # impose Hermiticity
        return d_matrix


# Monolayer dynamical matrix
class MonolayerDM(DM):
    def __init__(self, poscar_uc : Poscar, poscar_sc : Poscar, ph, GM_set, k_set, Gamma_idx, k_mags=None, using_flex=False):
        assert LA.norm(k_set[Gamma_idx]) == 0, f"Gamma index has nonzero norm {LA.norm(k_set[Gamma_idx])}"
        print("Symmetrizing force constants...")
        ph.symmetrize_force_constants()
        print("Force constants symmetrized.")
        self.n_at = sum(poscar_uc.natoms)
        self.uc = poscar_uc.structure; self.sc = poscar_sc.structure # get structure objects from Poscar objects
        self.GM_set = GM_set; self.n_GM = len(GM_set); self.k_set = k_set; self.ph = ph
        assert np.all(k_set[Gamma_idx] == [0,0]), f"Nonzero kGamma={k_set[Gamma_idx]}"
        self.n_k = len(k_set)
        self.pos_sc_id = self.__sc_atomic_id() if using_flex else None
        self.A0 = self.uc.lattice.matrix[:2,:2].T # remove z-axis
        self.G0 = 2 * pi * LA.inv(self.A0).T
        self.DM_set = None; self.dbgprint = True; self.l0_shape = None
        self.name = poscar_uc.comment; self.modes_built = False
        self.M = np.array([species.atomic_mass for species in poscar_uc.structure.species])
        if self.pos_sc_id is not None:
            print(f"MonolayerDM intralayer atomic IDs for {self.name}:", self.pos_sc_id)
        self.Gamma_idx = Gamma_idx
        self.k_mags = k_mags

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
        q = LA.inv(self.G0) @ q
        if len(q) != 3:
            q = np.append(q, np.zeros(3-len(q)))
        dm = ph.get_dynamical_matrix_at_q(q)
                    
        if self.dbgprint:
            print(f"Intralayer Level-0 shape: {dm.shape}")
            self.dbgprint = False
        if self.l0_shape is None:
            self.l0_shape = dm.shape
        return dm
    
    def plot_sampled_l0(self, outdir='.'):
        outdir = checkPath(os.path.abspath(outdir))
        lst = np.absolute([self.__block_intra_l0(GM, self.ph)[0,0] for GM in self.GM_set])
        plt.clf(); fig, ax = plt.subplots()
        ax.scatter(np.arange(self.n_GM), lst)
        fig.savefig(outdir + 'sampledl0.png')
        plt.close(fig)

    # Deprecated: Calculates dynamical matrix elements for direct or Cartesian coordinates `q`
    def __flex_block_intra_l0(self, q, ph):
        q = LA.inv(self.G0) @ q[:2]
        assert self.pos_sc_id is not None, "Fatal error in class initialization"
        pos_sc_idx = self.pos_sc_id
        smallest_vectors, multiplicity = ph.primitive.get_smallest_vectors()
        species = self.uc.species; uc_nat = len(species); d = 3 # Cartesian DOF
        ph.produce_force_constants()
        fc = ph.force_constants
        D = np.zeros([uc_nat*d, uc_nat*d], dtype=complex) # size: n_at x d (n_at of uc)

        """
        Iteration: choose a pair of atoms (given by sc ID, may be the same), fix one, 
        then iterate through all instances of the other. Add to the matrix the Fourier term
        divided by the multiplicity, which is determined by phonopy based on nearest-neighbor symmetry.
        """
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
                        vec = vec[0] * self.A0[:,0] + vec[1] * self.A0[:,1]
                        exp_here = np.exp(1j * np.dot(q, vec[0:2]))
                        d_matrix[id1:id1+3, id2:id2+3] += fc_here * exp_here / np.sqrt(m1) / np.sqrt(m2) / multi
        d_matrix = (d_matrix + d_matrix.conj().transpose()) / 2 # impose Hermiticity
        return d_matrix
    
    # Create level-1 block matrix (each matrix becomes a block in level-2)
    def __block_intra_l1(self):
        # print(self.k_set)
        # print("AFTER:")
        # kdir = self.k_set @ LA.inv(self.G0)
        # breakpoint()
        # ys = [LA.eigvals(self.__block_intra_l0(LA.inv(self.G0) @ k, self.ph)) for k in self.k_set]
        # plt.clf()
        # for k_mag, y in zip(self.k_mags, ys):
        #     plt.scatter([k_mag]*len(y), y, c='royalblue', s=0.07)
        # plt.xticks([])
        # plt.savefig("/Users/jonathanlu/Documents/tmos2/zoe/GOOSEBERRY.png")
        # sys.exit()
        self.__sc_atomic_id()

        self.DM_set = [block_diag(*[self.__flex_block_intra_l0(k+GM, self.ph) for GM in self.GM_set]) for k in self.k_set]
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
        return self.DM_set
    
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
        xlabs = (r'K', r'$\Gamma$', r'M')
        plt.xticks(corner_kmags, xlabs)
        plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
        title = "First-layer phonon modes"
        if name is not None:
            title += f" of {name}"
        plt.title(title)
        plt.savefig(outdir + filename)
        succ(f"Intralayer plot written to {outdir+filename}")
        return

    def get_GM_set(self):
        print(f"Retrieved GM sample set for solid {self.name}")
        return self.GM_set

    def get_k_set(self):
        print(f"Retrieved k sample set for solid {self.name}")
        return self.k_set


# Build interlayer dynamical matrix block via summing over configurations
class InterlayerDM(DM):
    def __init__(self, per_layer_at_idxs, bl_M, 
                 b_set, k_set, 
                 GM_set, G0_set, 
                 species_per_layer, 
                 ph_list=None, 
                 force_matrices=None):
        self.M = np.array([[species.atomic_mass for species in layer] for layer in species_per_layer])
        self.bl_M = bl_M; self.bl_n_at = len(bl_M)
        print(f"Interlayer DM uses bilayer with {self.bl_n_at} atoms of masses {self.bl_M}")
        def force_from_dm(i, dm):
            M = self.bl_M
            Mi = M[i]
            temp = np.split(dm[3*i:3*(i+1)], self.bl_n_at, axis=1)
            assert temp[0].shape == (3,3)
            return sum([sqrt(Mi*Mj) * subblk for subblk, Mj in zip(temp, M)])
        def force_adjust(f):
            for i in range(self.bl_n_at):
                blM = self.bl_M
                temp = np.split(f[3*i:3*(i+1)], self.bl_n_at, axis=1)
                temp = [sqrt(blM[j]/blM[i]) * blk for j, blk in enumerate(temp)]
                assert temp[0].shape == (3,3)
                f[3*i:3*(i+1), 3*i:3*(i+1)] -= sum(temp)
            return f

        assert (ph_list is not None) ^ (force_matrices is not None), "Must give exactly one of: phonopy obj list, force matrix list"
        assert len(b_set[0]) == 2, "Shift vectors must be 2-dimensional"
        self.b_set = b_set; self.ph_list = ph_list # list of phonopy objects for each config
        self.nshift = len(b_set); print(f"Number of configurations: {self.nshift}")
        assert int(sqrt(self.nshift))**2 == self.nshift, f"Number of shifts {self.nshift} must be a perfect square"
        self.GM_set = GM_set; self.G0_set = G0_set; self.DM = None
        self.per_layer_at_idxs = per_layer_at_idxs; assert len(self.per_layer_at_idxs) == 2, f"Only 2 layers supported"
        self.k_set = k_set
        self.modes_built = False
        self.n_GM = len(GM_set)
        if ph_list is not None:
            def sym_dm_at_gamma(i, ph):
                ph.symmetrize_force_constants()
                return ph.get_dynamical_matrix_at_q([0,0,0])
            self.force_matrices = [sym_dm_at_gamma(i, ph) for i, ph in enumerate(ph_list)] # DM(Gamma) ~= FC (mass-scaled)
        else:
            self.force_matrices = force_matrices
       
        assert self.force_matrices[0].shape[0] == self.force_matrices[0].shape[1], f"Force matrix is not square: shape {self.force_matrices[0].shape}"
        assert self.force_matrices[0].shape[0] % 2 == 0, f"Force matrix size is odd: shape {self.force_matrices[0].shape}"

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
        self.GMi_intra_blocks = [[0]*(n_GM-1) for i in range(2)]
        self.all_blocks = [0]*(n_GM)
        self.all_blocks[0] = D0
        self.Gamma_block = D0
        for i in range(1, n_GM): # fill first row/col
            self.DM[0][i], self.GMi_intra_blocks[0][i-1], self.GMi_intra_blocks[1][i-1] = self.__block_inter_l0(self.G0_set[i])
            self.all_blocks[i] = self.DM[0][i]
            self.DM[i][0],_,_ = self.__block_inter_l0(-self.G0_set[i])
            assert self.DM[0][i].shape == block_l0_shape and self.DM[i][0].shape == block_l0_shape, f"Shape GM0{i}={self.DM[0][i].shape}, GM{i}0={self.DM[i][0].shape}, expected {block_l0_shape}"
            assert np.isclose(LA.norm(self.DM[0][i]), LA.norm(self.DM[i][0]), rtol=1e-5), f"Level-0 interlayer DM blocks for G0{i} not inversion-symmetric:\n {LA.norm(self.DM[0][i])}\nvs. \n{LA.norm(self.DM[i][0])}"
        for i in range(n_GM): # fill diagonal
            self.DM[i][i] = D0
        self.DM = bmat(self.DM).toarray() # convert NoneTypes to zero-matrix blocks to make sparse matrix
        return self.DM
    
    def __intras_block_inter_l1(self):
        l1 = [0]*len(self.k_set); l2 = [0]*len(self.k_set)
        for i, k in enumerate(self.k_set):
            blocks = [self.__block_inter_l0(G0j, k=k) for G0j in self.G0_set]
            l1[i] = block_diag(*[b[1] for b in blocks])
            l2[i] = block_diag(*[b[2] for b in blocks])
        return l1, l2
    
    def get_intra_DM_set(self, k_mags=None, corner_kmags=None, outdir='.'):
        l1, l2 = self.__intras_block_inter_l1()
        if k_mags is not None and corner_kmags is not None:
            self.__plot_band(l1, l2, k_mags, corner_kmags, outdir=outdir)
        return l1, l2

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
    
    def __plot_band(self, l1, l2, k_mags, corner_kmags, outdir='./', 
                  filename='interintra_' + DEFAULT_PH_BAND_PLOT_NAME, name=None, cutoff=None):
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        if not self.modes_built:
            print("Configuration intralayer modes not built yet, building...")
            self.mode_set_per_layer = [self.__build_modes(k_mags, l, dump=False, outdir=outdir) for l in [l1, l2]]
            print("Configuration intralayer modes built.")
        plt.clf()
        for i, mode_set in enumerate(self.mode_set_per_layer):
            for k_mag, modes in mode_set:
                if cutoff is not None:
                    modes = modes[modes <= cutoff]
                plt.scatter([k_mag] * len(modes), modes, c='royalblue', s=0.07)
            xlabs = ('K', r'$\Gamma$', r'M')
            plt.xticks(corner_kmags, xlabs)
            plt.ylabel(r'$\omega\,(\mathrm{cm}^{-1})$')
            title = f"Layer-{i+1} phonon modes"
            if name is not None:
                title += f" of {name}"
            plt.title(title)
            this_filename = filename[:filename.index('.')] + f'_{i+1}' + filename[filename.index('.'):]
            plt.savefig(outdir + this_filename)
            print(f"Intralayer plot {i+1} written to {outdir+this_filename}")
    
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
        DMs_layer1 = l1.get_DM_set(); DMs_layer2 = l2.get_DM_set(); DM_inter = inter.get_DM()
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
        self.DMs = [self.__block_l2([DMs_layer1[i], DMs_layer2[i]], DM_inter) for i in range(self.n_k)]

    # Create level-2 (final level--full matrix) block matrix with intralayer and interlayer terms
    def __block_l2(self, DM_intras, DM_inter):
        assert len(DM_intras) == 2
        assert DM_intras[0].shape == DM_intras[1].shape == DM_inter.shape
        DM_intras = list(map(lambda D: (D + D.conjugate().T) / 2, DM_intras)) # impose Hermiticity
        DM = np.block([[DM_intras[0], DM_inter], [DM_inter.conjugate().T, DM_intras[1]]])
        return DM 
    
    # enforce an acoustic sum rule at Gamma
    def Gamma_acoustic_sum_rule(self):
        """
        Bugs: fails if number of atoms in layers differs. 
        Change `mid` to `self.szs[0]` to fix (easy).
        """
        update("ENFORCING: Gamma aoustic sum rule")
        def intra_adjust(l, i, j, l0sz, n_at, DMintra):
            assert l in [0,1]; Mi = self.M[l,i]
            idx = 3*i + j*l0sz
            temp = np.split(DMintra[idx:idx+3,j*l0sz:(j+1)*l0sz], n_at, axis=1)
            return sum([sqrt(Mj/Mi) * subblk for subblk, Mj in zip(temp, self.M[l])])
        def inter_adjust(l, i, n_at, blk): # n_at is number in *other* layer
            assert l in [0,1]; Mi = self.M[l,i]
            if l == 1:
                blk = blk.conjugate().T
            temp = np.split(blk[3*i:3*(i+1)], n_at, axis=1)
            return sum([sqrt(Mj/Mi) * subblk for subblk, Mj in zip(temp, self.M[1-l])])
        dm = self.get_DM_at_Gamma()
        l0szs = self.l0szs
        inter_Gamma = self.interobj.get_Gamma_block()
        for i in range(self.n_ats[0]): # layer 1 atoms
            intra = dm[:self.szs[0],:self.szs[0]]
            intrasum = np.real_if_close(intra_adjust(0, i, 0, l0szs[0], self.n_ats[0], intra))
            intersum = np.real_if_close(inter_adjust(0, i, self.n_ats[1], inter_Gamma))
            assert intrasum.shape == (3,3), f"Invalid intra shape {intrasum.shape}"
            assert intersum.shape == (3,3), f"Invalid inter shape {intersum.shape}"
            totalsum = np.real_if_close(intrasum + intersum)
            for k in range(self.n_k):
                self.DMs[k][3*i:3*(i+1),3*i:3*(i+1)] -= totalsum
            
        for i in range(self.n_ats[1]): # layer 2 atoms
            intra = dm[self.szs[0]:self.szs[0]+self.szs[1],self.szs[0]:self.szs[0]+self.szs[1]]
            assert dm.shape == tuple(map(lambda x: 2*x, intra.shape))
            intrasum = np.real_if_close(intra_adjust(1, i, 0, l0szs[1], self.n_ats[1], intra))
            intersum = np.real_if_close(inter_adjust(1, i, self.n_ats[0], inter_Gamma))
            assert intrasum.shape == (3,3), f"Invalid intra shape {intrasum.shape}"
            assert intersum.shape == (3,3), f"Invalid inter shape {intersum.shape}"
            totalsum = np.real_if_close(intrasum + intersum)
            mid = self.DMs[0].shape[0] // 2
            for j in range(self.n_k):
                self.DMs[j][mid + 3*i : mid + 3*(i+1), mid + 3*i : mid + 3*(i+1)] -= totalsum
    
    def Moire_acoustic_sum_rule(self, plotG=False):
        """
        Bugs: fails if number of atoms in layers differs. 
        Change `mid` to `self.szs[0]` to fix (easy).
        """
        update("ENFORCING: Moire off diagonal aoustic sum rule")
        def intra_adjust(l, i, j, l0sz, n_at, DMintra):
            assert l in [0,1]; Mi = self.M[l,i]
            idx = 3*i + j*l0sz
            temp = np.split(DMintra[idx:idx+3,j*l0sz:(j+1)*l0sz], n_at, axis=1)
            del temp[i] # remove diagonal II component
            M = [m for midx, m in enumerate(self.M[l]) if midx != i]
            return sum([sqrt(Mj/Mi) * subblk for subblk, Mj in zip(temp, M)])
        def inter_adjust(l, i, n_at, blk): # n_at is number in *other* layer
            assert l in [0,1]; Mi = self.M[l,i]
            if l == 1:
                blk = blk.conjugate().T
            temp = np.split(blk[3*i:3*(i+1)], n_at, axis=1)
            return sum([sqrt(Mj/Mi) * subblk for subblk, Mj in zip(temp, self.M[1-l])])
        def plot_Gpos(pltlist, i, cmpt=None, which='Intra', outdir='.'):
            outdir =checkPath( os.path.abspath(outdir))
            plt.clf(); fig, ax = plt.subplots()
            if cmpt is None:
                colors = list(map(lambda x: np.mean(np.absolute(x)), pltlist))
            else:
                colors = list(map(lambda x: x[cmpt], pltlist))
                cmptstr = ','.join(list(map(str, cmpt)))
            cf = ax.scatter(self.GM_set[1:,0], self.GM_set[1:,1], cmap='winter', c=colors)
            ax.set_xlabel(r'$\widetilde{k}_x$'); ax.set_ylabel(r'$\widetilde{k}_y$')
            plt.title(which + r"layer at-%d $\mathbf{\widetilde{G}} \neq 0$"%i + ("norm averages" if cmpt is None else f"for cmpt {cmpt}"))
            fig.colorbar(cf, ax=ax)
            outname = outdir + f'{which}_at{i}-{"avg" if cmpt is None else cmptstr}.png'
            fig.savefig(outname)
            print(f"Saved to {outname}")
            plt.close(fig)

        dm = self.get_DM_at_Gamma()
        l0szs = self.l0szs
        inter_Gammas = self.off_diag_blocks

        intra = dm[:self.szs[0],:self.szs[0]]
        for i in range(self.n_ats[0]): # layer 1 atoms
            intralist = [intra_adjust(0, i, j, l0szs[0], self.n_ats[0], intra) for j in range(1, self.n_GM)]
            interlist = [inter_adjust(0, i, self.n_ats[1], inter_Gamma) for inter_Gamma in inter_Gammas]
            intrasum = np.real_if_close(sum(intralist))
            intersum = np.real_if_close(sum(interlist))
            print(f"[At {i}] Intrasum:\n{intrasum}\nIntersum:\n{intersum}\n")
            assert intrasum.shape == (3,3), f"Invalid intra shape {intrasum.shape}"
            assert intersum.shape == (3,3), f"Invalid inter shape {intersum.shape}"
            totalsum = np.real_if_close(intrasum + intersum)
            print(f"At {i} at Gamma before Moire correction:\n{self.DMs[self.Gamma_idx][3*i:3*(i+1),3*i:3*(i+1)]}")
            for j in range(self.n_k):
                self.DMs[j][3*i:3*(i+1),3*i:3*(i+1)] -= intersum
            # Plot the matrix elements/averages
            print("PLOTTING: G != 0 components for intra and inter in layer 1 reference")
            if plotG:
                plot_Gpos(intralist, i, which='Intra')
                plot_Gpos(interlist, i, which='Inter')
                for j in [(1,1), (1,0), (0,2)]:
                    plot_Gpos(intralist, i, cmpt=j, which='Intra')
                    plot_Gpos(interlist, i, which='Inter', cmpt=j)
            print(f"At {i} at Gamma after Moire correction:\n{self.DMs[self.Gamma_idx][3*i:3*(i+1),3*i:3*(i+1)]}")
            
        intra = dm[self.szs[0]:self.szs[0]+self.szs[1],self.szs[0]:self.szs[0]+self.szs[1]]
        for i in range(self.n_ats[1]): # layer 2 atoms
            assert dm.shape == tuple(map(lambda x: 2*x, intra.shape))
            intrasum = np.real_if_close(sum([intra_adjust(1, i, j, l0szs[1], self.n_ats[1], intra) for j in range(1, self.n_GM)]))
            intersum = np.real_if_close(sum([inter_adjust(1, i, self.n_ats[0], inter_Gamma) for inter_Gamma in inter_Gammas]))
            assert intrasum.shape == (3,3), f"Invalid intra shape {intrasum.shape}"
            assert intersum.shape == (3,3), f"Invalid inter shape {intersum.shape}"
            totalsum = np.real_if_close(intrasum + intersum)
            mid = self.DMs[0].shape[0] // 2
            for j in range(self.n_k):
                self.DMs[j][mid + 3*i : mid + 3*(i+1), mid + 3*i : mid + 3*(i+1)] -= intersum
    
    # enforce a moire sum rule for off-diagonal blocks in interlayer component
    def dep_Moire_acoustic_sum_rule(self):
        # * On the first intralayer loop over the 3x3 block diagonal i and sum over all 3x3 blocks
        # * in the interlayer 1-2 matrix in row i. Repeat for second intralayer and interlayer 2-1.
        update("ENFORCING: deprecated moire acoustic sum rule")
        print(f"Moire acoustic sum rule on: {len(self.off_diag_blocks)} GM vectors")
        M = self.M
        mid = self.DMs[0].shape[0] // 2; print(f"Slicing intra at midpoint index {mid}")
        for k in range(self.n_k):
            for i in range(self.n_ats[0]): 
                Mi = M[0,i]
                # intralayer 1
                self.DMs[k][3*i:3*(i+1), 3*i:3*(i+1)] -= sum([sum([sqrt(M[0,j//3] / Mi) * Gblk[3*i:3*(i+1),j:j+3] for j in range(0,Gblk.shape[1],3) if j != 3*i]) for Gblk in self.on_diag_blocks[0]])
                # interlayer 1-2
                self.DMs[k][3*i:3*(i+1), 3*i:3*(i+1)] -= sum([sum([sqrt(M[1,j//3] / Mi) * Gblk[3*i:3*(i+1),j:j+3] for j in range(0,Gblk.shape[1],3)]) for Gblk in self.off_diag_blocks])
            for i in range(self.n_ats[1]): # interlayer 2-1
                Mi = M[1,i]
                # intralayer 2
                self.DMs[k][mid+3*i : mid+3*(i+1), mid+3*i : mid+3*(i+1)] -= sum([sum([sqrt(M[1,j//3] / Mi) * Gblk[3*i:3*(i+1),j:j+3] for j in range(0,Gblk.shape[1],3) if j != 3*i]) for Gblk in self.on_diag_blocks[1]])
                # interlayer 2-1
                self.DMs[k][mid+3*i : mid+3*(i+1), mid+3*i : mid+3*(i+1)] -= sum([sum([sqrt(M[0,j//3] / Mi) * Gblk.conjugate().T[3*i:3*(i+1),j:j+3] for j in range(0,Gblk.shape[0],3)]) for Gblk in self.off_diag_blocks])

    # Retrieve list dynamical matrices corresponding to the list of sampled k-vectors
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


