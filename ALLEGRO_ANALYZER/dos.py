"""
Frontend script for computing the phonon dispersion, 
realspace transforms, and Gamma modes in 
configuration space for a twisted cell.
"""
import numpy as np
import numpy.linalg as LA
from math import sqrt
from time import time
from phonopy import Phonopy
import phonopy
import pymatgen.core.structure as struct
from pymatgen.io.vasp.inputs import Poscar
from __class_DM import TwistedDM, InterlayerDM, MonolayerDM, TwistedDOS
from __class_PhonopyAPI import PhonopyAPI
from bzsampler import get_bz_sample
from ___constants_sampling import DEFAULT_KDIM
from ___constants_names import (
    SPOSCAR_NAME, PH_FORCE_SETS_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME, 
    POSCAR_CONFIG_NAMEPRE, 
    SHIFT_NAME, SHIFTS_NPY_NAME, 
    DEFAULT_PH_BAND_PLOT_NAME, 
    ANGLE_SAMPLE_INAME, 
    MODE_TNSR_ONAME, ANGLE_SAMPLE_ONAME, GAMMA_IDX_ONAME, 
    K_MAGS_ONAME, K_SET_ONAME, DM_TNSR_ONAME
)
from ___constants_compute import DEFAULT_NUM_LAYERS
from __directory_searchers import checkPath, findDirsinDir, findFilesInDir
from __dirModifications import build_dir
from __class_Configuration import Configuration
from __class_RelaxerAPI import RelaxerAPI
from __class_ForceInterp import ForceInterp, FourierForceInterp
from __class_PhononConfig import TwistedRealspacePhonon
from ___helpers_parsing import greet, update, succ, warn, err, is_flag, check_not_flag
import os, sys, copy
from math import pi
    
if __name__ == '__main__':
    start = time()
    USAGE_MSG = f"Usage: python3 {sys.argv[0]} -deg <twist angle (deg)> -name <solid name> -cut <frequency cutoff> -dir <base I/O dir> -o <output dir>"
    args = sys.argv[1:]; i = 0; n = len(args)
    indir = '.'; outdir = '.'; theta = None; name = None; outname = DEFAULT_PH_BAND_PLOT_NAME; width = None
    plot_intra = False; relax = False; realspace = False; force_sum = False; multirelax = False
    do_sum_rule = True; kdim = DEFAULT_KDIM
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
            print(f"theta: {np.rad2deg(theta)}")
        elif args[i] == '--ns':
            i += 1; do_sum_rule = False; print("NOT using sum rule...")
        elif args[i] == '-name':
            i += 1; check_not_flag(args[i]); name = args[i]; i += 1
        elif args[i] == '-k':
            i += 1; check_not_flag(args[i]); kdim = int(args[i])
            assert kdim > 20, f"Mesh size {kdim} too small, try at least 21"
            i += 1
        elif args[i] == '--width':
            i += 1; check_not_flag(args[i]); width = float(args[i]); i += 1
        elif args[i] == '-fout':
            i += 1; check_not_flag(args[i]); outname = args[i]; i += 1
        elif args[i] == '--intra':
            i += 1; plot_intra = True
        elif args[i] == '--relax' or args[i] == '-r':
            i += 1
            relax = True
        elif args[i] == '--mr':
            i += 1
            multirelax = True
        elif args[i] in ['--usage', '--help', '-h']:
            print(USAGE_MSG)
            sys.exit(0)
        else:
            warn(f'Warning: unknown flag "{args[i]}" ignored')
            i += 1
    if not theta:
        err(f"Error: must supply twist angle. Run `python3 {sys.argv[0]} --usage` for help.")
    indir = checkPath(os.path.abspath(indir)); outdir = checkPath(os.path.abspath(outdir))
    assert os.path.isdir(indir), f"Directory {indir} does not exist"
    assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
    if multirelax:
        assert not relax, f"`multirelax` and `relax` are incompatible options"
        assert os.path.isfile(indir + ANGLE_SAMPLE_INAME), f"Must provide input file {ANGLE_SAMPLE_INAME}"
    print(f"WD: {indir}, Output to: {outdir}")
    outname = 'd' + str(round(np.rad2deg(theta))) + "_" + outname
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
    n_at = [len(posc.structure.species) for posc in poscars_uc]
    A0 = s0.T
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
    print(f"Found monolayer supercell dimensions to be {super_dim}")

    update("Retrieving per-layer atomic indices...")
    per_layer_at_idxs = Configuration.load_at_idxs(build_dir([indir, CONFIG_DIR_NAME]))
    assert len(per_layer_at_idxs) == 2, f"Only 2 layers supported (for now), got {len(per_layer_at_idxs)}"
    print(f"Found per-layer atomic indices.")
    
    update("Generating phonopy objects via API...")
    ph_api = PhonopyAPI(indir) # retreive phonopy objects
    ml_ph_list = ph_api.intra_ph_list()
    print("Phonopy objects generated.")

    update("Sampling G and k sets...")
    bzsamples = get_bz_sample(theta, poscars_uc[0], outdir, make_plot=True, super_dim=super_dim)
    _, GM_set = bzsamples.get_GM_set(); _, G0_set = bzsamples.get_G0_set()
    k_mesh = bzsamples.sample_mesh_k(kdim=kdim)
    bzsamples.plot_mesh_sampling()
    print("Sampling complete.")
    os.chdir(outdir)

    update("Constructing intralayer dynamical matrix objects...")
    MLDMs = [MonolayerDM(uc, sc, ph, GM_set, G0_set, k_mesh, 0) \
                         for uc, sc, ph in zip(poscars_uc, poscars_sc, ml_ph_list)]
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
    b_set = np.array(b_set)
    b_set = b_set @ A0.T # to cartesian
    print("Shift vectors imported.")
        
    print("Getting interlayer phonopy objects from API...")
    config_ph_list = ph_api.inter_ph_list()
    print("Phonopy objects retrieved.")

    if multirelax:
        thetas = None
        with open(indir + ANGLE_SAMPLE_INAME, 'r') as f:
            thetas = list(map(float, f.read().splitlines()))
            assert len(thetas) == 3 and 20 > thetas[1] > thetas[0] > 0 and thetas[2] > 0, f"Invalid thetas {thetas}"
        ntheta = int(thetas[2])
        thetas = np.linspace(np.deg2rad(thetas[0]), np.deg2rad(thetas[1]), ntheta)
        # TODO DOS on each theta

    else:
        relaxed_forces = None; b_relaxed = None
        if relax:
            print("Non-uniformizing configurations via relaxation...")
            relax_api = RelaxerAPI(np.rad2deg(theta), gridsz, outdir, s0.T)
            b_relaxed = relax_api.get_configs(cartesian=True)
            relax_api.plot_relaxation()

        print("Note: Using GM sampling set from intralayer calculations.")
        print("Constructing interlayer dynamical matrix objects...")
        bl_M = ph_api.bl_masses()
        print(f"Bilayer masses: {list(bl_M)}")
        ILDM = InterlayerDM(per_layer_at_idxs, bl_M, 
                             b_set, 
                             bzsamples.get_kpts0()[0], 
                             GM_set, G0_set, 
                             [p.structure.species for p in poscars_uc], 
                             ph_list=config_ph_list, 
                             force_matrices=None, br_set=b_relaxed)
        print("Interlayer DM objects constructed.")

        print("Combining into a single twisted dynamical matrix object...")
        TDM = TwistedDM(MLDMs[0], MLDMs[1], ILDM, np.zeros(len(k_mesh)), 
                        [p.structure.species for p in poscars_uc], 0)
        if do_sum_rule:
            TDM.apply_sum_rule()
        print("Twisted dynamical matrix object constructed.")

        TDMs = TDM.get_DM_set()
        TDMs_intra = TDM.get_intra_set()
        dos = TwistedDOS(TDMs, TDMs_intra, len(GM_set), np.rad2deg(theta), width=width, kdim=kdim)
        print("DOS object constructed")
        dos.plot_DOS(vertical=False, outdir=outdir)

    succ("Successfully completed phonon DOS analysis (Took %.3lfs)."%(time()-start))


