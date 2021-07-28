# Dynamical matrix calculations for twisted materials
import numpy as np
import sys, copy, os
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
from pymatgen.io.vasp.inputs import Poscar
import numpy.linalg as LA
from math import pi, floor
from itertools import product as prod
from time import time
import matplotlib.pyplot as plt
from ___constants_names import DIR_SHORTCUTS
from ___constants_sampling import *
from __class_BZSampler import BZSampler
from ___helpers_parsing import greet, succ, warn, err, is_flag, check_not_flag
from __directory_searchers import checkPath
from ___constants_phonopy import SUPER_DIM
from ____exit_with_error import exit_with_error

def get_bz_sample(theta, poscar, outdir, super_dim=SUPER_DIM[:2], 
                  max_shell=DEFAULT_MAX_SHELL, gridsz=DEFAULT_GRID_SIZE, nk=DEFAULT_NK, 
                  log=True, make_plot=False):
    assert isinstance(poscar, Poscar) or isinstance(poscar, str), "Input 'poscar' must either be a 'Poscar' obj or a path to a POSCAR file"
    s = poscar.structure if isinstance(poscar, Poscar) else struct.Structure.from_file(poscar)
    if isinstance(poscar, str):
        assert os.path.isfile(poscar), f"Invalid path {poscar}"
    assert os.path.isdir(outdir), f"Invalid directory {outdir}"
    assert max_shell >= 1 and gridsz > 1 and nk > 1 and (theta is None or theta > 0), "Invalid sampling parameters"

    # Build realspace A-basis and reciprocal space G-basis in the moire cell
    A0 = s.lattice.matrix[0:2, 0:2].T # makes realspace lattice A-basis matrix [a1 a2]
    G0 = 2 * pi * LA.inv(A0).T 
    if log:
        print("A0:\n", A0)
        print("G0:\n", G0)

    sample = BZSampler(A0, G0, theta, outdir=outdir, log=log, super_dim=super_dim)
    sample.sample_GM(mn_grid_sz=gridsz, max_shell=max_shell)
    sample.sample_k(nk=nk, log=True)
    sample.sample_k0(nk=nk, log=True)
    if make_plot:
        sample.plot_sampling()
        sample.plot_sampling0()
    return sample

if __name__ == '__main__':
    # Parse input
    start = time()
    USAGE_MSG = '\nUsage: python3 <DIR>/bzsampler.py -deg <twist degree> -dir <main dir> -o <output dir> -f <POSCAR name>\n\n\t(Optional flags: -gsz <G grid size> -nk <kpts per line> -sh <max shell>)\n'
    args = sys.argv[1:]; i = 0; n = len(args)
    theta = None; indir = '.'; outdir = None; pname = 'POSCAR'; p_found = False; sample_G0 = True
    gridsz = DEFAULT_GRID_SIZE; max_shell = DEFAULT_MAX_SHELL; nk = DEFAULT_NK
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-deg':
            i += 1; check_not_flag(args[i]); theta = np.deg2rad(float(args[i])); i += 1
            sample_G0 = False
        elif args[i] == '-dir':
            i += 1; check_not_flag(args[i])
            indir = checkPath(os.path.abspath(args[i]))
            if args[i] in DIR_SHORTCUTS or args[i][0] == '.':
                warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
            i += 1
        elif args[i] == '-o':
            i += 1; check_not_flag(args[i])
            outdir = checkPath(os.path.abspath(args[i]))
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

    if not indir:
        err(USAGE_MSG)
    elif theta is not None and 180 < theta < 0:
        err(f"Error: invalid twist angle {theta} degrees")
    elif not outdir:
        outdir = indir
        warn(f"No output directory specified. Defaulting to input directory {indir}...")
    if not p_found:
        warn("No POSCAR name given as input. Defaulting to name 'POSCAR'...")
    assert os.path.isdir(indir), f'Invalid input directory {indir}'
    assert os.path.isdir(outdir), f'Invalid output directory {outdir}'
    assert os.path.isfile(indir+pname), f'No POSCAR with name {pname} found in {indir}'
    greet("DM calculator starting...")
    print(f"Sampling {'pristine' if sample_G0 else 'moire'} reciprocal lattice vectors...")
    if not sample_G0:
        print("Twist angle: {np.rad2deg(theta)} deg/{'%.3lf'%theta} rad")
    print(f"WD: {indir}\ngrid size = {gridsz}, max shell = {max_shell}, nk = {nk}")

    get_bz_sample(theta, indir + pname, outdir, 
                  max_shell=max_shell, gridsz=gridsz, nk=nk, 
                  log=True, make_plot=True)
                  
    succ("Sampling successfully completed (took %.3lfs)."%(time()-start))

