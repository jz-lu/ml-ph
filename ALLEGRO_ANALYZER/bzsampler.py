# Dynamical matrix calculations for twisted materials
import numpy as np
import sys, copy, os
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
import numpy.linalg as LA
from math import pi, floor
from itertools import product as prod
from time import time, sleep
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
        R1 = np.array([[np.cos(theta/2), -np.sin(theta/2)], [np.sin(theta/2), np.cos(theta/2)]])
        R2 = LA.inv(R1)
        G1 = np.matmul(R1, G0); G2 = np.matmul(R2, G0) 
        self.GM = G1 - G2 # moire G-basis
        self.AM = LA.inv(self.GM).T / (2 * pi) # moire A-basis
        if log:
            print("Moire G-basis:\n", self.GM)
            print("Moire A-basis:\n", self.AM)
        self.g_idxs = None; self.g_arr = None; self.kline = None; self.kmags = None; self.corners = None
        self.ltype = lattice_type
        if lattice_type != 'hexagonal':
            exit_with_error("Error: k-points sampler does not (yet) support non-hexagonal lattices")

    # Sample Gtilde vectors
    def sample_G(self, mn_grid_sz=DEFAULT_GRID_SIZE, max_shell=DEFAULT_MAX_SHELL):
        grid = np.arange(-mn_grid_sz, mn_grid_sz + 1); self.nsh = max_shell
        GM1 = self.GM[:,0]; GM2 = self.GM[:,1] # Moire reciprocal lattice vectors
        g_arr = np.array([m*GM1 + n*GM2 for m, n in prod(grid, grid)])
        g_idxs = np.array(list(prod(grid, grid)))

        # Filter out any G that does not satisfy the closeness condition |GM| < (shell+0.1) * |GM1|
        g_cutoff = (floor(max_shell)+0.1) * LA.norm(GM1) # add 0.1 for numerical error
        cut_makers = (np.sqrt(g_arr[:,0]**2 + g_arr[:,1]**2) <= g_cutoff) # indicators for which G made the cutoff
        g_idxs = g_idxs[cut_makers,:]
        g_arr = g_arr[cut_makers,:] # special numpy syntax for filtering by an indicator array `cut_makers`
        self.g_idxs = g_idxs; self.g_arr = g_arr
        return (g_idxs, g_arr)

    # Sample nk points along each IBZ boundary line
    def sample_k(self, nk=DEFAULT_NK, log=False):
        assert self.ltype == 'hexagonal'
        G00 = G0[:,0]; G01 = G0[:,1]; d = G00.shape[0]
        Gamma = np.zeros(d); K = 1/3 * (G00 + G01); M = 1/2 * G00
        ncorners = 4; corners = np.zeros([ncorners, d]) # k-points IBZ boundary corners
        corners[0,:] = Gamma; corners[1,:] = K; corners[2,:] = M; corners[3,:] = Gamma
        self.corners = corners

        # Sample nk-1 (drop last point) per line (ncorners-1 lines)
        nsample = (nk-1) * (ncorners-1)
        kline = np.zeros([nsample, d]); kmags = np.zeros(nsample); kmag_start = LA.norm(Gamma)
        corner_kmags = [kmag_start]
        for line in range(ncorners-1): # skip last point equals first, so skip it
            kidx = line*(nk-1) # convert line index to k-index
            # Drop second corner point in each line to avoid resampling corner points
            kline[kidx : kidx+nk-1] = np.linspace(corners[line], corners[line+1], nk)[:-1]
            dline_mag = LA.norm(corners[line+1] - corners[line])
            mags = np.linspace(kmag_start, kmag_start + dline_mag, nk)
            kmag_start = mags[-1] # update start point of magnitude to end of current line
            corner_kmags.append(kmag_start)
            kmags[kidx : kidx+nk-1] = mags[:-1]
        self.kline = kline; self.kmags = kmags
        if log:
            print("Corner magnitudes:", corner_kmags)
        return (kline, kmags)

    def plot_sampling(self, filename='sampling.png'):
        assert self.g_idxs is not None, "Cannot plot sampling until G vectors have been sampled"
        assert self.corners is not None, "Cannot plot sampling until k-points have been sampled"
        labels = (r'$\Gamma$', r'K', r'M', r'$\Gamma$')
        plt.clf()
        _, ax = plt.subplots()
        plt.scatter(self.g_arr[:,0], self.g_arr[:,1], 
                    c='black', 
                    label=r'$\widetilde{\mathbf{G}}_{mn}$ in shell %d'%self.nsh)
        plt.plot(self.kline[:,0], self.kline[:,1], 
                 c='teal', 
                 label=r'$\mathbf{k}$-points (%d pts)'%len(self.kline))
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

if __name__ == '__main__':
    # Parse input
    USAGE_MSG = '\nUsage: python3 <DIR>/bzsampler.py -deg <twist degree> -dir <main dir> -o <output dir> -f <POSCAR name>\n\n\t(Optional flags: -gsz <G grid size> -nk <kpts per line> -sh <max shell>)\n'
    args = sys.argv[1:]; i = 0; n = len(args)
    theta = None; indir = '.'; outdir = None; pname = 'POSCAR'; p_found = False
    gridsz = DEFAULT_GRID_SIZE; max_shell = DEFAULT_MAX_SHELL; nk = DEFAULT_NK
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-deg':
            i += 1; check_not_flag(args[i]); theta = np.deg2rad(float(args[i])); i += 1
        elif args[i] == '-dir':
            i += 1; check_not_flag(args[i])
            indir = checkPath(args[i])
            if args[i] in DIR_SHORTCUTS or args[i][0] == '.':
                warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
            i += 1
        elif args[i] == '-o':
            i += 1; check_not_flag(args[i])
            outdir = checkPath(args[i])
            if args[i] in DIR_SHORTCUTS or args[i][0] == '.':
                warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
            i += 1
        elif args[i] in ['-f', '-p']:
            assert not p_found, 'Multiple POSCAR filenames given'
            p_found = True
            i += 1; check_not_flag(args[i])
            pname = args[i]
            i += 1
        elif args[i] == '-gsz':
            i += 1; check_not_flag(args[i]); gridsz = int(args[i]); i += 1
        elif args[i] == '-sh':
            i += 1; check_not_flag(args[i]); max_shell = int(args[i]); i += 1
        elif args[i] == '-nk':
            i += 1; check_not_flag(args[i]); nk = int(args[i]); i += 1
        else:
            err(f"Error: unknown token {args[i]}")

    if not (theta and indir):
        err(USAGE_MSG)
    elif 180 < theta < 0:
        err(f"Error: invalid twist angle {theta} degrees")
    elif not outdir:
        outdir = indir
        warn(f"No output directory specified. Defaulting to input directory {indir}...")
    if not p_found:
        warn("No POSCAR name given as input. Defaulting to name 'POSCAR'...")
    assert os.path.isdir(indir), f'Invalid input directory {indir}'
    assert os.path.isdir(outdir), f'Invalid output directory {indir}'
    assert os.path.isfile(indir + pname), f'No POSCAR with name {pname} found in {indir}'
    greet("DM calculator starting...")
    print(f"Twist angle: {np.rad2deg(theta)} deg/{'%.3lf'%theta} rad\nWD: {indir}\ngrid size = {gridsz}, max shell = {max_shell}, nk = {nk}")

    # Build realspace A-basis and reciprocal space G-basis in the moire cell
    s = struct.Structure.from_file(indir + pname)
    A0 = s.lattice.matrix[0:2, 0:2].T # makes realspace lattice A-basis matrix [a1 a2]
    G0 = 2 * pi * LA.inv(A0).T 
    print("A0:\n", A0)
    print("G0:\n", G0)

    sample = BZSampler(G0, theta, outdir=outdir, log=True)
    sample.sample_G(mn_grid_sz=gridsz, max_shell=max_shell)
    sample.sample_k(nk=nk, log=True)
    sample.plot_sampling()
    succ("Sampling successfully completed.")

