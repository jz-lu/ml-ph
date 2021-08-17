# Sample moire G space and k-line space
import numpy as np
import sys, copy, os
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
import numpy.linalg as LA
from math import pi, floor
from itertools import product as prod
from time import time
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from ____exit_with_error import exit_with_error
from ___constants_names import DIR_SHORTCUTS
from ___constants_sampling import *
from ___constants_phonopy import SUPER_DIM
from ___helpers_parsing import greet, succ, warn, err, is_flag, check_not_flag
from __class_CarCollector import CarCollector
from __directory_searchers import checkPath
import random

class BZSampler:
    def __init__(self, A0, G0, theta=None, outdir='./', log=False, lattice_type='hexagonal', super_dim=SUPER_DIM[:2]):
        G0 = np.array(G0)
        theta = theta if theta is not None else 1
        assert pi > theta >= 0 # in rad
        assert os.path.isdir(outdir)
        self.outdir = checkPath(outdir); self.theta = theta; self.G0 = G0; self.A0 = A0
        self.lattice_angle = CarCollector.lattice_basis_angle(self.A0, col=True)
        self.super_dim = super_dim

        # Generate G basis via rotation matrices
        R1 = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]]) # rot by theta/2
        R2 = LA.inv(R1) # rot by -theta/2
        G1 = np.matmul(R1, G0); G2 = np.matmul(R2, G0); self.G0 = G0
        self.GM = G1 - G2 # moire G-basis
        self.AM = LA.inv(self.GM).T * 2 * pi # moire A-basis
        if log:
            print(f"Lattice angle: {self.lattice_angle}")
            print("Moire G-basis:\n", self.GM)
            print("Moire A-basis:\n", self.AM)
        self.g_idxs = None; self.GM_set = None; self.G0_set = None
        self.k_set = None; self.kmags = None; self.corners = None
        self.ltype = lattice_type
        self.GM_sampled = False; self.k_sampled = False; self.k0_sampled = False
        self.nsh = None
        if lattice_type != 'hexagonal':
            exit_with_error("Error: k-points sampler does not (yet) support non-hexagonal lattices")

    # Sample Gtilde vectors
    def sample_GM(self, mn_grid_sz=DEFAULT_GRID_SIZE, max_shell=DEFAULT_MAX_SHELL, tol=0.1):
        assert self.theta is not None, f"Must provide twist angle to sample moire G vectors"
        grid = np.arange(-mn_grid_sz, mn_grid_sz + 1); self.nsh = max_shell
        G01 = self.G0[:,0]; G02 = self.G0[:,1] # untwisted monolayer reciprocal lattice vectors
        GM1 = self.GM[:,0]; GM2 = self.GM[:,1] # Moire reciprocal lattice vectors
        GM_set = np.array([m*GM1 + n*GM2 for m, n in prod(grid, grid)])
        G0_set = np.array([m*G01 + n*G02 for m, n in prod(grid, grid)])
        g_idxs = np.array(list(prod(grid, grid)))
        self.nsh = max_shell; self.tol = tol

        # Filter out any G that does not satisfy the closeness condition |GM| < (shell+0.1) * |GM1|
        g_cutoff = (floor(max_shell) + tol) * LA.norm(GM1) # add 0.1 for numerical error
        cut_makers = (np.sqrt(GM_set[:,0]**2 + GM_set[:,1]**2) <= g_cutoff) # indicators for which G made the cutoff
        self.g_idxs = g_idxs[cut_makers,:]; self.GM_set = GM_set[cut_makers,:]; self.G0_set = G0_set[cut_makers,:]
        normsort = lambda x: LA.norm(x[0])
        self.g_idxs = np.array([g for _,g in sorted(zip(self.GM_set, self.g_idxs), key=normsort)]) # sort by norm of GM
        self.G0_set = np.array([g for _,g in sorted(zip(self.GM_set, self.G0_set), key=normsort)])
        self.GM_set = np.array(sorted(self.GM_set, key=LA.norm))
        print(f"GM set: \n{self.GM_set}\nG0 set: \n{self.G0_set}\nG indices: \n{self.g_idxs}")
        self.GM_sampled = True
        return (g_idxs, GM_set)

    # Sample Gtilde vectors
    def sample_G0(self, mn_grid_sz=DEFAULT_GRID_SIZE, max_shell=DEFAULT_MAX_SHELL, tol=0.1):
        grid = np.arange(-mn_grid_sz, mn_grid_sz + 1); self.nsh = max_shell
        G01 = self.G0[:,0]; G02 = self.G0[:,1] # untwisted monolayer reciprocal lattice vectors
        G0_set = np.array([m*G01 + n*G02 for m, n in prod(grid, grid)])
        g_idxs = np.array(list(prod(grid, grid)))

        # Filter out any G that does not satisfy the closeness condition |GM| < (shell+0.1) * |GM1|
        g_cutoff = (floor(max_shell) + tol) * LA.norm(G01) # add 0.1 for numerical error
        cut_makers = (np.sqrt(G0_set[:,0]**2 + G0_set[:,1]**2) <= g_cutoff) # indicators for which G made the cutoff
        g_idxs = g_idxs[cut_makers,:]; G0_set = G0_set[cut_makers,:]
        normsort = lambda x: LA.norm(x[0])
        g_idxs = np.array([g for _,g in sorted(zip(G0_set, g_idxs), key=normsort)]) # sort by norm of GM
        G0_set = np.array(sorted(self.GM_set, key=LA.norm))
        print(f"G0 set: \n{self.G0_set}\nG indices: \n{self.g_idxs}")
        return (g_idxs, G0_set)

    # Sample nk points along each IBZ boundary line
    def sample_k(self, nk=DEFAULT_NK, log=False):
        assert self.ltype == 'hexagonal'
        GM1 = self.GM[:,0]; GM2 = self.GM[:,1] # Moire reciprocal lattice vectors
        d = GM1.shape[0]
        print("k-sampler using K-Gamma-M-K sequence defined by moire reciprocal lattice vectors GM")
        Gamma = np.zeros(d); K = 1/3 * (GM1 + GM2); M = 1/2 * GM1
        if np.isclose(self.lattice_angle, 60):
            print("Using 60-degree unit cell for BZ...")
            M += 1/2 * GM2
            K += 1/3 * GM1
        else:
            print("Using 120-degree unit cell for BZ...")
        
        assert np.isclose(self.lattice_angle, 60) or np.isclose(self.lattice_angle, 120), f"k-sampler expects lattice angle to be 120 or 60 deg, but is {self.lattice_angle}"
        ncorners = 4; corners = np.zeros([ncorners, d]) # k-points IBZ boundary corners
        corners[0,:] = K; corners[1,:] = Gamma; corners[2,:] = M; corners[3,:] = K
        self.corners = corners

        # Sample nk-1 (drop last point) per line (ncorners-1 lines)
        nsample = (nk-1) * (ncorners-1)
        k_set = np.zeros([nsample, d]); kmags = np.zeros(nsample); kmag_start = 0
        corner_kmags = []
        for line in range(ncorners-1): # last point equals first, so skip it
            kidx = line*(nk-1) # convert line index to k-index
            # Drop second corner point in each line to avoid resampling corner points
            k_set[kidx : kidx+nk-1] = np.linspace(corners[line], corners[line+1], nk)[:-1]
            dline_mag = LA.norm(corners[line+1] - corners[line])
            mags = np.linspace(kmag_start, kmag_start + dline_mag, nk)
            corner_kmags.append(kmag_start)
            kmag_start = mags[-1] # update start point of flattened-k to end of current line
            kmags[kidx : kidx+nk-1] = mags[:-1]
        self.Gamma_idx = nk-1
        self.k_set = k_set; self.kmags = kmags
        self.k_sampled = True
        if log:
            print("Corner magnitudes:", corner_kmags)
            print(f"k at Gamma index {self.Gamma_idx}: {k_set[self.Gamma_idx]}")
        self.corner_kmags = corner_kmags
        return (k_set, kmags)
    
    # Sample nk points along each IBZ boundary line
    def sample_k0(self, nk=DEFAULT_NK, log=False):
        assert self.ltype == 'hexagonal'
        G01 = self.G0[:,0]; G02 = self.G0[:,1] # Pristine reciprocal lattice vectors
        d = G01.shape[0]
        print("k-sampler using Gamma-K-M sequence defined by pristine reciprocal lattice vectors G0")
        Gamma = np.zeros(d); K = 1/3 * (G01 + G02); M = 1/2 * G01
        if np.isclose(self.lattice_angle, 60):
            print("Using 60-degree unit cell for BZ...")
            M += 1/2 * G02
            K += 1/3 * G01
        else:
            print("Using 120-degree unit cell for BZ...")
        
        assert np.isclose(self.lattice_angle, 60) or np.isclose(self.lattice_angle, 120), f"k-sampler expects lattice angle to be 120 or 60 deg, but is {self.lattice_angle}"
        ncorners = 4; corners = np.zeros([ncorners, d]) # k-points IBZ boundary corners
        corners[0,:] = K; corners[1,:] = Gamma; corners[2,:] = M; corners[3,:] = K
        self.corners0 = corners

        # Sample nk-1 (drop last point) per line (ncorners-1 lines)
        nsample = (nk-1) * (ncorners-1)
        k_set = np.zeros([nsample, d]); kmags = np.zeros(nsample); kmag_start = 0
        corner_kmags = []
        for line in range(ncorners-1): # last point equals first, so skip it
            kidx = line*(nk-1) # convert line index to k-index
            # Drop second corner point in each line to avoid resampling corner points
            k_set[kidx : kidx+nk-1] = np.linspace(corners[line], corners[line+1], nk)[:-1]
            dline_mag = LA.norm(corners[line+1] - corners[line])
            mags = np.linspace(kmag_start, kmag_start + dline_mag, nk)
            corner_kmags.append(kmag_start)
            kmag_start = mags[-1] # update start point of flattened-k to end of current line
            kmags[kidx : kidx+nk-1] = mags[:-1]
        self.k0_set = k_set; self.k0mags = kmags
        self.corner_kmags0 = corner_kmags
        self.k0_sampled = True
        return (k_set, kmags)
    
    def plot_sampling0(self, filename='sampling0.png'):
        assert self.g_idxs is not None, "Cannot plot sampling until G vectors have been sampled"
        assert self.corners0 is not None, "Cannot plot sampling until k-points have been sampled"
        labels = (r'K', r'$\Gamma$', r'M', r'K')
        plt.clf()
        _, ax = plt.subplots()
        cols = list(mcolors.TABLEAU_COLORS.keys()); ncol = len(cols)
        random.shuffle(cols)
        G01 = self.G0_set[:,0]; G02 = self.G0_set[:,1]
        last_cut = 0
        for i in range(1, self.nsh+1):
            g_cutoff = (i + self.tol) * LA.norm(self.G0[:,0])
            cut_makers = np.logical_and(np.sqrt(G01**2 + G02**2) >= last_cut, np.sqrt(G01**2 + G02**2) <= g_cutoff)
            last_cut = g_cutoff
            plt.scatter(G01[cut_makers], G02[cut_makers], 
                        c=cols[i%ncol], 
                        label=r'$\widetilde{\mathbf{G}}_{mn}$ in shell %d'%i)
        for idx, G0 in zip(self.g_idxs, self.G0_set):
            plt.annotate(str(tuple(idx)), G0, # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(np.sign(G0[0])*10, -np.sign(G0[1])*10), # distance from text to points (x,y)
                 ha='center', fontsize=9.5) # horizontal alignment can be left, right or center
        plt.plot(self.k0_set[:,0], self.k0_set[:,1], 
                 c='teal', 
                 label=r'$\mathbf{k}$-points (%d pts)'%len(self.k_set))
        ax.set_aspect('equal') # prevent stretching of space in plot
        for (x, y), lab in zip(self.corners0, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-10,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        plt.title("Pristine Sampling Space")
        ax.set_xlabel(r'$k_x$'); ax.set_ylabel(r'$k_y$')
        outname = self.outdir + filename
        plt.savefig(outname); succ("Successfully wrote sampling plot out to " + outname)

    def plot_sampling(self, filename='sampling.png'):
        assert self.g_idxs is not None, "Cannot plot sampling until G vectors have been sampled"
        assert self.corners is not None, "Cannot plot sampling until k-points have been sampled"
        labels = (r'K', r'$\Gamma$', r'M', r'K')
        plt.clf()
        _, ax = plt.subplots()
        cols = list(mcolors.TABLEAU_COLORS.keys()); ncol = len(cols)
        random.shuffle(cols)
        GM1 = self.GM_set[:,0]; GM2 = self.GM_set[:,1]
        last_cut = 0
        for i in range(1, self.nsh+1):
            g_cutoff = (i + self.tol) * LA.norm(self.GM[:,0])
            cut_makers = np.logical_and(np.sqrt(GM1**2 + GM2**2) >= last_cut, np.sqrt(GM1**2 + GM2**2) <= g_cutoff)
            last_cut = g_cutoff
            plt.scatter(GM1[cut_makers], GM2[cut_makers], 
                        c=cols[i%ncol], 
                        label=r'$\widetilde{\mathbf{G}}_{mn}$ in shell %d'%i)
        for idx, GM in zip(self.g_idxs, self.GM_set):
            plt.annotate(str(tuple(idx)), GM, # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(np.sign(GM[0])*10, -np.sign(GM[1])*10), # distance from text to points (x,y)
                 ha='center', fontsize=9.5) # horizontal alignment can be left, right or center
        plt.plot(self.k_set[:,0], self.k_set[:,1], 
                 c='teal', 
                 label=r'$\mathbf{k}$-points (%d pts)'%len(self.k_set))
        ax.set_aspect('equal') # prevent stretching of space in plot
        for (x, y), lab in zip(self.corners, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-10,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        # plt.legend()
        plt.title("Sampling Space")
        ax.set_xlabel(r'$k_x$'); ax.set_ylabel(r'$k_y$')
        outname = self.outdir + filename
        plt.savefig(outname); succ("Successfully wrote sampling plot out to " + outname)

    def get_GM_set(self):
        assert self.GM_sampled, "Must run GM-sampler before retrieving GM-set"
        return self.g_idxs, self.GM_set
    
    def get_G0_set(self):
        assert self.GM_sampled, "Must run GM-sampler before retrieving G0-set"
        return self.g_idxs, self.G0_set

    def get_kpts(self):
        assert self.k_sampled, "Must run k-sampler before retrieving k-set"
        return self.k_set, self.kmags
    
    def get_kpts0(self):
        assert self.k0_sampled, "Must run k-sampler before retrieving k-set"
        return self.k0_set, self.k0mags
    
    def get_Gamma_idx(self):
        assert self.k_sampled, "Must run k-sampler before retrieving k-set"
        return self.Gamma_idx

    def get_corner_kmags(self):
        return self.corner_kmags

    def get_corner_kmags0(self):
        return self.corner_kmags0
    
    def sampled(self):
        return self.GM_sampled and self.k_sampled


