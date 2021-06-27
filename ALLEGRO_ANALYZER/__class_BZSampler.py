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
from ___constants_names import DIR_SHORTCUTS
from ___constants_sampling import *
from ___helpers_parsing import greet, succ, warn, err, is_flag, check_not_flag
from __directory_searchers import checkPath
from ____exit_with_error import exit_with_error

class BZSampler:
    def __init__(self, G0, theta, outdir='./', log=False, lattice_type='hexagonal'):
        G0 = np.array(G0) # ensure proper typing
        assert 180 > theta > 0
        assert os.path.isdir(outdir)
        self.outdir = checkPath(outdir); self.theta = theta; self.G0 = G0

        # Generate G basis via rotation matrices
        R1 = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]]) # rot by theta/2
        R2 = LA.inv(R1) # rot by -theta/2
        G1 = np.matmul(R1, G0); G2 = np.matmul(R2, G0); self.G0 = G0
        self.GM = G1 - G2 # moire G-basis
        self.AM = LA.inv(self.GM).T / (2 * pi) # moire A-basis
        if log:
            print("Moire G-basis:\n", self.GM)
            print("Moire A-basis:\n", self.AM)
        self.g_idxs = None; self.GM_set = None; self.G0_set = None
        self.k_set = None; self.kmags = None; self.corners = None
        self.ltype = lattice_type
        self.GM_sampled = False; self.k_sampled = False
        if lattice_type != 'hexagonal':
            exit_with_error("Error: k-points sampler does not (yet) support non-hexagonal lattices")

    # Sample Gtilde vectors
    def sample_GM(self, mn_grid_sz=DEFAULT_GRID_SIZE, max_shell=DEFAULT_MAX_SHELL, tol=0.1):
        grid = np.arange(-mn_grid_sz, mn_grid_sz + 1); self.nsh = max_shell
        G01 = self.G0[:,0]; G02 = self.G0[:,1] # untwisted monolayer reciprocal lattice vectors
        GM1 = self.GM[:,0]; GM2 = self.GM[:,1] # Moire reciprocal lattice vectors
        GM_set = np.array([m*GM1 + n*GM2 for m, n in prod(grid, grid)])
        G0_set = np.array([m*G01 + n*G02 for m, n in prod(grid, grid)])
        g_idxs = np.array(list(prod(grid, grid)))

        # Filter out any G that does not satisfy the closeness condition |GM| < (shell+0.1) * |GM1|
        g_cutoff = (floor(max_shell) + tol) * LA.norm(GM1) # add 0.1 for numerical error
        cut_makers = (np.sqrt(GM_set[:,0]**2 + GM_set[:,1]**2) <= g_cutoff) # indicators for which G made the cutoff
        self.g_idxs = g_idxs[cut_makers,:]
        self.GM_set = GM_set[cut_makers,:] # special numpy syntax for filtering by an indicator array `cut_makers`
        self.G0_set = G0_set[cut_makers,:]
        print(f"GM set: \n{self.GM_set}\nG0 set: \n{self.G0_set}\nG indices: \n{self.g_idxs}")
        self.GM_sampled = True
        return (g_idxs, GM_set)

    # Sample nk points along each IBZ boundary line
    def sample_k(self, nk=DEFAULT_NK, log=False):
        assert self.ltype == 'hexagonal'
        G00 = self.G0[:,0]; G01 = self.G0[:,1]; d = G00.shape[0]
        Gamma = np.zeros(d); K = 1/3 * (G00 + G01); M = 1/2 * G00
        ncorners = 4; corners = np.zeros([ncorners, d]) # k-points IBZ boundary corners
        corners[0,:] = Gamma; corners[1,:] = K; corners[2,:] = M; corners[3,:] = Gamma
        self.corners = corners

        # Sample nk-1 (drop last point) per line (ncorners-1 lines)
        nsample = (nk-1) * (ncorners-1)
        k_set = np.zeros([nsample, d]); kmags = np.zeros(nsample); kmag_start = LA.norm(Gamma)
        corner_kmags = []
        for line in range(ncorners-1): # last point equals first, so skip it
            kidx = line*(nk-1) # convert line index to k-index
            # Drop second corner point in each line to avoid resampling corner points
            k_set[kidx : kidx+nk-1] = np.linspace(corners[line], corners[line+1], nk)[:-1]
            dline_mag = LA.norm(corners[line+1] - corners[line])
            mags = np.linspace(kmag_start, kmag_start + dline_mag, nk)
            corner_kmags.append(kmag_start)
            kmag_start = mags[-1] # update start point of magnitude to end of current line
            kmags[kidx : kidx+nk-1] = mags[:-1]
        corner_kmags += [corner_kmags[0]]
        self.k_set = k_set; self.kmags = kmags
        self.k_sampled = True
        if log:
            print("Corner magnitudes:", corner_kmags)
        self.corner_kmags = corner_kmags
        return (k_set, kmags)

    def plot_sampling(self, filename='sampling.png'):
        assert self.g_idxs is not None, "Cannot plot sampling until G vectors have been sampled"
        assert self.corners is not None, "Cannot plot sampling until k-points have been sampled"
        labels = (r'$\Gamma$', r'K', r'M', r'$\Gamma$')
        plt.clf()
        _, ax = plt.subplots()
        plt.scatter(self.GM_set[:,0], self.GM_set[:,1], 
                    c='black', 
                    label=r'$\widetilde{\mathbf{G}}_{mn}$ in shell %d'%self.nsh)
        plt.plot(self.k_set[:,0], self.k_set[:,1], 
                 c='teal', 
                 label=r'$\mathbf{k}$-points (%d pts)'%len(self.k_set))
        ax.set_aspect('equal') # prevent stretching of space in plot
        for (x, y), lab in zip(self.corners, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-10,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        plt.legend(); plt.title("Sampling Space")
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

    def get_corner_kmags(self):
        return self.corner_kmags
    
    def sampled(self):
        return self.GM_sampled and self.k_sampled


