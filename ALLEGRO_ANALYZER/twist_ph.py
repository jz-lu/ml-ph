# Compute the phonon modes for a twisted cell
import numpy as np
import numpy.linalg as LA
from math import sqrt
from time import time
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
from pymatgen.io.vasp.inputs import Poscar
from __class_DM import TwistedDM, InterlayerDM, MonolayerDM
from __class_PhonopyAPI import PhonopyAPI
from bzsampler import get_bz_sample
from ___constants_names import (
    SPOSCAR_NAME, PH_FORCE_SETS_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME, 
    POSCAR_CONFIG_NAMEPRE, 
    SHIFT_NAME, SHIFTS_NPY_NAME, 
    DEFAULT_PH_BAND_PLOT_NAME
)
from ___constants_compute import DEFAULT_NUM_LAYERS
from __directory_searchers import checkPath, findDirsinDir, findFilesInDir
from __dirModifications import build_dir
from __class_Configuration import Configuration
from __class_RelaxerAPI import RelaxerAPI
from __class_ForceInterp import ForceInterp, FourierForceInterp
from ___helpers_parsing import greet, update, succ, warn, err, is_flag, check_not_flag
import os, sys

RELAX_FOURIER_INTERP = 1
RELAX_SPLINE_INTERP = 2

if __name__ == '__main__':
    start = time()
    USAGE_MSG = f"Usage: python3 {sys.argv[0]} -deg <twist angle (deg)> -name <solid name> -cut <frequency cutoff> -dir <base I/O dir> -o <output dir>"
    args = sys.argv[1:]; i = 0; n = len(args)
    indir = '.'; outdir = '.'; theta = None; name = None; outname = DEFAULT_PH_BAND_PLOT_NAME; cutoff = None
    plot_intra = False; relax = False
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-dir':
            i += 1; check_not_flag(args[i]); indir = args[i]; i += 1
        elif args[i] == '-o':
            i += 1; check_not_flag(args[i]); outdir = args[i]; i += 1
        elif args[i] == '-deg':
            i += 1; check_not_flag(args[i]); theta = np.deg2rad(float(args[i])); i += 1
        elif args[i] == '-name':
            i += 1; check_not_flag(args[i]); name = args[i]; i += 1
        elif args[i] == '-cut':
            i += 1; check_not_flag(args[i]); cutoff = float(args[i]); i += 1
        elif args[i] == '-fout':
            i += 1; check_not_flag(args[i]); outname = args[i]; i += 1
        elif args[i] == '--intra':
            i += 1; plot_intra = True
        elif args[i] == '--relax' or args[i] == 'r':
            i += 1
            relax = RELAX_SPLINE_INTERP if args[i] in ['r', 'real', 'R', 'REAL'] else RELAX_FOURIER_INTERP
            i += 1 
        elif args[i] in ['--usage', '--help']:
            print(USAGE_MSG)
            sys.exit(0)
        else:
            warn(f'Warning: unknown flag "{args[i]} ignored')
            i += 1
    if not theta:
        err(f"Error: must supply twist angle. Run `python3 {sys.argv[0]} --usage` for help.")
    indir = checkPath(os.path.abspath(indir)); outdir = checkPath(os.path.abspath(outdir))
    assert os.path.isdir(indir), f"Directory {indir} does not exist"
    assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
    print(f"WD: {indir}, Output to: {outdir}")
    outname = 'd' + str(int(np.rad2deg(theta))) + "_" + outname
    print("Building twisted crystal phonon modes...")
    warn("Note: you must have already ran twisted calculation in given directory and retrieved forces from phonopy")

    greet("Working on intralayer components...")
    print("Searching for POSCAR inputs...")
    poscars_uc = sorted(findFilesInDir(indir, POSCAR_CONFIG_NAMEPRE, searchType='start')); npos = len(poscars_uc)
    poscars_uc = [Poscar.from_file(indir + pname) for pname in poscars_uc]
    assert npos > 1, "Must give at least 2 POSCARs, 1 per layer, for twist calculations"
    assert npos == DEFAULT_NUM_LAYERS, "Twist calculations for more than 2 layers not supported (yet)"
    s0 = poscars_uc[0].structure.lattice.matrix[0:2, 0:2]
    s1 = poscars_uc[1].structure.lattice.matrix[0:2, 0:2]
    assert np.allclose(s0, s1), "Input POSCARs must have the same lattice matrices"
    print("POSCAR inputs found.")

    print("Searching for SPOSCAR inputs...")
    poscars_sc = [''] * npos
    for i in range(npos):
        sposcar_path = build_dir([indir, MONOLAYER_DIR_NAME + str(i+1), ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME]) + SPOSCAR_NAME
        assert os.path.isfile(sposcar_path), f"SPOSCAR not found at {sposcar_path}"
        poscars_sc[i] = Poscar.from_file(sposcar_path)
    print("SPOSCAR inputs found.")
    
    S0 = poscars_sc[0].structure.lattice.matrix[:2,:2]
    Xcol = S0[:,0]; Ycol = S0[:,1]; xcol = s0[:,0]; ycol = s0[:,1]
    xscale = [round(X/x) for x, X in zip(xcol, Xcol) if not np.isclose(x, 0)][0]
    yscale = [round(Y/y) for y, Y in zip(ycol, Ycol) if not np.isclose(y, 0)][0]
    super_dim = (xscale, yscale)
    print(f"Found supercell dimensions to be {super_dim}")

    update("Retrieving per-layer atomic indices...")
    per_layer_at_idxs = Configuration.load_at_idxs(build_dir([indir, CONFIG_DIR_NAME]))
    assert len(per_layer_at_idxs) == 2, f"Only 2 layers supported (for now), got {len(per_layer_at_idxs)}"
    print(f"Found per-layer atomic indices: {per_layer_at_idxs}")
    
    update("Generating phonopy objects via API...")
    ph_api = PhonopyAPI(indir) # retreive phonopy objects
    ml_ph_list = ph_api.intra_ph_list()
    print("Phonopy objects generated.")

    update("Sampling G and k sets...")
    bzsamples = get_bz_sample(theta, poscars_uc[0], outdir, make_plot=True, super_dim=super_dim)
    _, GM_set = bzsamples.get_GM_set(); _, G0_set = bzsamples.get_G0_set()
    k_set, k_mags = bzsamples.get_kpts(); corner_kmags = bzsamples.get_corner_kmags()
    print("Sampling complete.")

    update("Constructing intralayer dynamical matrix objects...")
    MLDMs = [MonolayerDM(uc, sc, ph, GM_set, k_set) for uc, sc, ph in zip(poscars_uc, poscars_sc, ml_ph_list)]
    assert len(MLDMs) == 2, f"Only 2-layer solids supported for twisted DM (for now), got {len(MLDMs)}"
    if plot_intra:
        print("Plotting one intralayer component...")
        MLDMs[0].plot_band(k_mags, corner_kmags, name=name, outdir=outdir, filename='intra_'+outname, cutoff=cutoff)
    print("Intralayer DM objects constructed.")

    greet("Working on interlayer components...")
    print("Importing shift vectors...")
    nshift = len(findDirsinDir(build_dir([indir, CONFIG_DIR_NAME]), CONFIG_SUBDIR_NAME, searchType='start'))
    gridsz = int(sqrt(nshift)); assert gridsz**2 == nshift, f"Number of shifts {nshift} must be a perfect square"
    b_set = ['']*nshift
    for i in range(nshift):
        shift_file_path = build_dir([indir, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME + str(i)]) + SHIFT_NAME
        assert os.path.isfile(shift_file_path), f"{shift_file_path} not found"
        with open(shift_file_path) as f:
            b_set[i] = np.array(list(map(float, f.read().splitlines()))[:2]) # shifts must be 2-dimensional
    print("Shift vectors imported:", b_set)
        

    print("Getting interlayer phonopy objects from API...")
    config_ph_list = ph_api.inter_ph_list()
    print("Phonopy objects retrieved.")

    relaxed_forces = None
    if relax:
        print("Non-uniformizing configurations via relaxation...")
        relax_api = RelaxerAPI(np.rad2deg(theta), gridsz, outdir, s0.T)
        b_relaxed = relax_api.get_configs()
        if relax == RELAX_SPLINE_INTERP:
            print("Using spline interpolation...")
            interp = ForceInterp(config_ph_list, np.array(b_set))
            relaxed_forces = interp.fc_tnsr_at(b_relaxed)
        else:
            print("Using Fourier interpolation...")
            interp = FourierForceInterp(config_ph_list, np.array(b_set), s0.T)
            relaxed_forces = interp.predict(b_relaxed)
        print(f"Relaxed interpolation force tensor is of shape: {relaxed_forces.shape}")
        # Relaxed forces tensor has index (f1, f2, b), but needs to be of form (b, f1, f2)
        relaxed_forces = np.transpose(relaxed_forces, axes=(2, 0, 1))
        print(f"Transposed to shape: {relaxed_forces.shape}")

    print("Note: Using GM sampling set from intralayer calculations.")
    print("Constructing interlayer dynamical matrix objects...")
    ILDM = None
    config_ph_list = None if relax else config_ph_list
    ILDM =  InterlayerDM(per_layer_at_idxs, b_set, GM_set, G0_set, ph_list=config_ph_list, force_matrices=relaxed_forces)
    evals = LA.eigvals(ILDM.get_DM())
    signs = (-1)*(evals < 0) + (evals > 0) # pull negative sign out of square root to plot imaginary frequencies
    modes_k = signs * np.sqrt(np.abs(evals)) * (15.633302*33.356) # eV/Angs^2 -> THz ~ 15.633302; THz -> cm^-1 ~ 33.356
    print("Interlayer modes (k-independent):", modes_k[modes_k != 0])
    print("Interlayer DM objects constructed.")

    print("Combining into a single twisted dynamical matrix object...")
    TDM = TwistedDM(MLDMs[0], MLDMs[1], ILDM, k_mags)
    print("Twisted dynamical matrix object constructed.")

    print(f"Diagonalizing and outputting modes with corners {corner_kmags}...")
    TDM.plot_band(corner_kmags, np.rad2deg(theta), outdir=outdir, name=name, filename=outname, cutoff=cutoff)
    print("Modes outputted.")

    succ("Successfully completed phonon mode analysis (Took %.3lfs)."%(time()-start))


