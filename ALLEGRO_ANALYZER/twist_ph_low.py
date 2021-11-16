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
from __class_DM import TwistedDM, InterlayerDM, MonolayerDM, TwistedDOS, TwistedPlotter
from __class_PhonopyAPI import PhonopyAPI
from bzsampler import get_bz_sample
from ___constants_names import (
    SPOSCAR_NAME, PH_FORCE_SETS_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME, 
    POSCAR_CONFIG_NAMEPRE, 
    SHIFT_NAME, SHIFTS_NPY_NAME, 
    DEFAULT_PH_BAND_PLOT_NAME, DEFAULT_PH_BANDDOS_PLOT_NAME, 
    ANGLE_SAMPLE_INAME, RSPC_LIST_INAME, RSPC_K_INAME, 
    MODE_TNSR_ONAME, ANGLE_SAMPLE_ONAME, GAMMA_IDX_ONAME, 
    K_MAGS_ONAME, K_SET_ONAME, DM_TNSR_ONAME
)
from ___constants_phonopy import POSCAR_UNIT_NAME
from ___constants_compute import DEFAULT_NUM_LAYERS
from ___constants_sampling import IBZ_GAMMA, IBZ_K_60, IBZ_M_60, DEFAULT_RSPC_SZ
from __directory_searchers import checkPath, findDirsinDir, findFilesInDir
from __dirModifications import build_dir
from __class_Configuration import Configuration
from __class_RelaxerAPI import RelaxerAPI
from __class_ForceInterp import ForceInterp, FourierForceInterp
from __class_PhononConfig import TwistedRealspacePhonon
from ___helpers_parsing import greet, update, succ, warn, err, is_flag, check_not_flag
import os, sys, copy, argparse
from math import pi, sqrt

#* deprecated
RELAX_FOURIER_INTERP = 1
RELAX_SPLINE_INTERP = 2 
    
if __name__ == '__main__':
    start = time()
    parser = argparse.ArgumentParser(description="Moire phonon analysis machine")
    parser.add_argument("theta", type=float, help='twist angle')
    parser.add_argument("-d", "--dir", type=str, help='main directory path', default=".")
    parser.add_argument("-o", "--out", type=str, help='output path', default=".")
    parser.add_argument("--ns", action="store_true", help='do not apply sum rule')
    parser.add_argument("-n", "--name", type=str, help='material name')
    parser.add_argument("-c", "--cut", type=int, help='band range cutoff')
    parser.add_argument("--oname", type=str, help='output filename', default=DEFAULT_PH_BANDDOS_PLOT_NAME)
    parser.add_argument("--intra", action="store_true", help='do intralayer analysis')
    parser.add_argument("-r", "--relax", action="store_true", help='do theta-relaxation')
    parser.add_argument("--mr", action="store_true", help=f'multirelax (req parameter file \'{ANGLE_SAMPLE_INAME}\')')
    parser.add_argument("--rs", action="store_true", help='do realspace analysis')
    parser.add_argument("--kdir", action="store_true", help="input k-points are in direct coordinates (default Cartesian)")
    parser.add_argument("--rssz", type=int, help="supercell size for realspace plots", default=DEFAULT_RSPC_SZ)
    parser.add_argument("--dos", type=int, help='DOS k-mesh size')
    parser.add_argument("-f", "--fsum", action="store_true", help='output force sums (only useful for debugging)')
    args = parser.parse_args()
    
    # Legacy compatibility
    theta = np.deg2rad(args.theta/10); indir = args.dir; outdir = args.out; multirelax = args.mr
    relax = args.relax; plot_intra = args.intra; force_sum = args.fsum; name = args.name
    cutoff = args.cut; realspace = args.rs; do_sum_rule = not args.ns; outname = args.oname
    print(f"Twist angle: {round(np.rad2deg(theta), 6)} deg")

    if cutoff is None:
        cutoff = int(40 * ((args.theta/10)**0.85))

    if not theta:
        err(f"Error: must supply twist angle. Run `python3 {sys.argv[0]} --usage` for help.")
    indir = checkPath(os.path.abspath(indir)); outdir = checkPath(os.path.abspath(outdir))
    assert os.path.isdir(indir), f"Directory {indir} does not exist"
    assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
    if multirelax:
        assert not relax, f"`multirelax` and `relax` are incompatible options"
        assert os.path.isfile(indir + ANGLE_SAMPLE_INAME), f"Must provide input file {ANGLE_SAMPLE_INAME}"
    if realspace:
        assert os.path.isfile(indir + RSPC_LIST_INAME), f"Must provide input file {RSPC_LIST_INAME}"
        assert args.rssz > 0, f"Realspace supercell size must be a positive integer, but got {args.rssz}"
    print(f"WD: {indir}, Output to: {outdir}")
    outname = 'd' + str(round(np.rad2deg(theta), 6)) + "_" + outname
    print("Building twisted crystal phonon modes...")
    warn("Note: you must have already ran twisted calculation in given directory and retrieved forces from phonopy")

    greet("Working on intralayer components...")
    print("Searching for POSCAR inputs...")
    poscars_uc = sorted(findFilesInDir(indir, POSCAR_CONFIG_NAMEPRE, searchType='start')); npos = len(poscars_uc)
    poscars_uc = [Poscar.from_file(indir + pname) for pname in poscars_uc]
    assert npos > 1, "Must give at least 2 POSCARs, 1 per layer, for twist calculations"
    assert npos == DEFAULT_NUM_LAYERS, "Twist calculations for more than 2 layers not supported (yet)"
    assert (args.dos is None) or (args.dos > 0), f"Must give a positive k-mesh size"
    if (args.dos is not None) and (args.dos < 21):
        warn(f"WARNING: realspace k-mesh size {args.dos} likely too small to be accurate")
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

    print("Searching for config (S)POSCAR inputs...")
    cfgpath = build_dir([indir, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME+"0", \
        ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
    cfg_sc = Poscar.from_file(cfgpath + SPOSCAR_NAME).structure
    cfg_uc = Poscar.from_file(cfgpath + POSCAR_UNIT_NAME).structure
    print("Config (S)POSCAR inputs found.")
    
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
    Gamma_idx = bzsamples.get_Gamma_idx()
    k_set, k_mags = bzsamples.get_kpts(); corner_kmags = bzsamples.get_corner_kmags()
    k0_set, k0_mags = bzsamples.get_kpts0(); corner_k0mags = bzsamples.get_corner_kmags0()
    # np.save(outdir + "km.npy", k_mags)
    # np.save(outdir + "k.npy", k_set)
    # np.save(outdir + "ck.npy", corner_kmags)
    print("Sampling complete.")
    os.chdir(outdir)

    update("Constructing intralayer dynamical matrix objects...")
    MLDMs = [MonolayerDM(uc, sc, ph, GM_set, \
                         G0_set, k_set, Gamma_idx, k_mags=k_mags) \
                         for uc, sc, ph in zip(poscars_uc, poscars_sc, ml_ph_list)]
    if plot_intra:
        print("Plotting one intralayer component...")
        print("Plotting in pristine space...")
        MLDMs[0].plot_pristine_band(k0_set, k0_mags, corner_k0mags, outdir=outdir)
        print("Plotting in moire space...") # ; MLDMs[0].plot_sampled_l0()
        MLDMs[0].plot_band(k_mags, corner_kmags, name=name, outdir=outdir, filename='intra_'+outname, cutoff=cutoff)
    if force_sum:
        for i, MLDM in enumerate(MLDMs):
            print(f"Computing layer {i} force sums...")
            MLDM.print_force_sum()
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
    b_set = b_set @ A0.T
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
        dmat_dim = 3 * len(GM_set) * sum(n_at); print(f"Moire dynamical matrix size = {dmat_dim}")
        dmat_tnsr = np.zeros((ntheta, len(k_set), dmat_dim, dmat_dim))
        # TODO wrong rn!! Need to resample k for each theta!!
        for i, theta in enumerate(thetas):
            # See non-multirelax case for comments
            relax_api = RelaxerAPI(round(np.rad2deg(theta), 6), gridsz, outdir, s0.T)
            b_relaxed = relax_api.get_configs(cartesian=True)

            interp = FourierForceInterp(config_ph_list, np.array(b_set), s0.T)
            relaxed_forces = interp.predict(b_relaxed)
            relaxed_forces = np.transpose(relaxed_forces, axes=(4,0,1,2,3))

            relaxed_ph_list = copy.deepcopy(config_ph_list)
            for i in range(len(config_ph_list)):
                relaxed_ph_list[i].force_constants = relaxed_forces[i]

            bl_M = ph_api.bl_masses()
            if i == 0:
                print(f"Bilayer masses: {bl_M}")
            ILDM =  InterlayerDM(per_layer_at_idxs, bl_M, 
                                 b_relaxed, k_set, 
                                 GM_set, G0_set, 
                                 [p.structure.species for p in poscars_uc], cfg_sc, cfg_uc, 
                                 ph_list=relaxed_ph_list, 
                                 force_matrices=None)
            TDM = TwistedDM(MLDMs[0], MLDMs[1], ILDM, k_mags, [p.structure.species for p in poscars_uc], Gamma_idx)
            TDM.apply_sum_rule()
            dmat_tnsr[i] = TDM.get_DM_set()
            # np.save(outdir + MODE_TNSR_ONAME%i, TDM.k_mode_tensor())
        np.save(outdir + ANGLE_SAMPLE_ONAME, thetas)
        np.save(outdir + GAMMA_IDX_ONAME, Gamma_idx)
        np.save(outdir + K_MAGS_ONAME, k_mags)
        np.save(outdir + K_SET_ONAME, k_set)
        np.save(outdir + DM_TNSR_ONAME, dmat_tnsr)
        update(f"Saved all multirelax output files to {outdir}")

    else:
        relaxed_forces = None; b_relaxed = None
        if relax:
            print("Non-uniformizing configurations via relaxation...")
            relax_api = RelaxerAPI(round(np.rad2deg(theta), 6), gridsz, outdir, s0.T, dedensify=True) # keep dense sampling for fourier interp
            b_relaxed = relax_api.get_configs(cartesian=True)
            print(f"UPDATED: configuration grid size to {int(sqrt(b_relaxed.shape[0]))} for implicit Fourier interpolation")
            relax_api.plot_relaxation()

        print("Note: Using GM sampling set from intralayer calculations.")
        print("Constructing interlayer dynamical matrix objects...")
        bl_M = ph_api.bl_masses()
        print(f"Bilayer masses: {list(bl_M)}")
        ILDM = InterlayerDM(per_layer_at_idxs, bl_M, 
                             b_set, 
                             bzsamples.get_kpts0()[0], 
                             GM_set, G0_set, 
                             [p.structure.species for p in poscars_uc], cfg_sc, cfg_uc, 
                             ph_list=config_ph_list, 
                             force_matrices=None, br_set=b_relaxed, A0=A0)
        print("Interlayer DM objects constructed.")

        print("Combining into a single twisted dynamical matrix object...")
        TDM = TwistedDM(MLDMs[0], MLDMs[1], ILDM, k_mags, [p.structure.species for p in poscars_uc], Gamma_idx)
        print("Twisted dynamical matrix object constructed.") 
        if do_sum_rule:
            TDM.apply_sum_rule()
        mode_set = TDM.get_mode_set()
        TDM.plot_band(corner_kmags, round(np.rad2deg(theta), 6), outdir=outdir, name=name, \
            filename='band_'+outname, cutoff=cutoff)

        if force_sum:
            print("\nAnalyzing sum of forces in twisted bilayer (G0 only)...")
            TDM.print_force_sum()

        if realspace: 
            print("Starting realspace analysis...")
            rspc_kpts = np.array([IBZ_GAMMA]) # default to Gamma point
            if os.path.isfile(indir + RSPC_K_INAME):
                rspc_kpts = np.concatenate((rspc_kpts, np.loadtxt(indir + RSPC_K_INAME)))
                print(f"Loaded list of k-points from {indir + RSPC_K_INAME}")
            else:
                print("No input set of k-points found, using Gamma point as default")
            rspc_kpts = np.unique(rspc_kpts, axis=0) # filter duplicate rows
            if args.kdir:
                print("USING: k-points in direct coordinates")
                dir_rspc_kpts = rspc_kpts
                rspc_kpts = bzsamples.k_dir_to_cart(rspc_kpts) # convert to Cartesian coordinates
            else:
                print("USING: k-points in Cartesian coordinates")
                dir_rspc_kpts = bzsamples.k_cart_to_dir(rspc_kpts) # convert to direct coordinates
            num_kpts = rspc_kpts.shape[0]
            assert rspc_kpts.shape[1] == 2
            print(f"Analyzing phonons in realspace at k-points ::\nCartesian:\n{rspc_kpts}\nDirect:\n{dir_rspc_kpts}")
            n_at = sum([sum(p.natoms) for p in poscars_uc])
            print(f"Number of atoms in bilayer configuration cell: {n_at}")

            # Build the dynamical matrices at the requested k-points
            rspc_MLDMs = [MonolayerDM(uc, sc, ph, GM_set, G0_set, rspc_kpts, 0) \
                        for uc, sc, ph in zip(poscars_uc, poscars_sc, ml_ph_list)]
            rspc_TDM = TwistedDM(rspc_MLDMs[0], rspc_MLDMs[1], ILDM, np.zeros(rspc_kpts.shape[0]), \
                [p.structure.species for p in poscars_uc], 0)
            if do_sum_rule:
                rspc_TDM.apply_sum_rule()
            rspc_TDMs = rspc_TDM.get_DM_set()
            modeidxs = np.loadtxt(indir + RSPC_LIST_INAME)
            for kidx, rspc_k in enumerate(rspc_kpts):
                dir_rspc_k = dir_rspc_kpts[kidx]
                kpt_name = "(%4.3lf, %4.3lf)"%(dir_rspc_k[0], dir_rspc_k[1])
                if np.isclose(IBZ_GAMMA, dir_rspc_k).all():
                    kpt_name = r'$\Gamma$'
                elif np.isclose(IBZ_K_60, dir_rspc_k).all():
                    kpt_name = 'K'
                elif np.isclose(IBZ_M_60, dir_rspc_k).all():
                    kpt_name = 'M'
                log_name = kpt_name[2:-1] if kpt_name == r'$\Gamma$' else "%4.3lf_%4.3lf"%(dir_rspc_k[0], dir_rspc_k[1])
                this_outdir = outdir
                if num_kpts > 1:
                    this_outdir += log_name
                    os.mkdir(this_outdir)
                print(f"[{kidx+1}/{num_kpts}] NOW WORKING ON: k = {log_name}, print to {this_outdir}")
                twrph = TwistedRealspacePhonon(round(np.rad2deg(theta), 6), rspc_k, GM_set, 
                        rspc_TDMs[kidx], n_at, bl_M, poscars_uc, outdir=this_outdir, modeidxs=modeidxs, 
                        kpt=kpt_name, RSPC_SUPERCELL_SIZE=args.rssz)
                twrph.plot_phonons()
                print(f"Phonons in realspace analyzed at {rspc_k} [Cart], {dir_rspc_k} [dir] (i.e. {kpt_name}).")
                # twrph.plot_spatial_avgs()
                # twrph.plot_atomic_avgs()
        
        print(f"Diagonalizing and outputting modes with corners {corner_kmags}...")
        if args.dos is not None:
            print("Including DOS...")
            k_mesh = bzsamples.sample_mesh_k(kdim=args.dos)
            bzsamples.plot_mesh_sampling()
            mesh_MLDMs = [MonolayerDM(uc, sc, ph, GM_set, G0_set, k_mesh, 0) \
                         for uc, sc, ph in zip(poscars_uc, poscars_sc, ml_ph_list)]
            mesh_TDM = TwistedDM(mesh_MLDMs[0], mesh_MLDMs[1], ILDM, np.zeros(len(k_mesh)), 
                        [p.structure.species for p in poscars_uc], 0)
            if do_sum_rule:
                mesh_TDM.apply_sum_rule()
            mesh_TDMs = mesh_TDM.get_DM_set()
            mesh_TDMs_intra = mesh_TDM.get_intra_set()
            widths = np.linspace(0.01, 0.21, 21)
            if np.rad2deg(theta) > 1.9:
                widths = np.linspace(0.11, 0.51, 21)
            eigsys = None
            for WIDTH in widths:
            # WIDTH = 0.05
                # iDOS = TwistedDOS(mesh_TDMs_intra, len(GM_set), np.rad2deg(theta), width=WIDTH, kdim=args.dos)
                # normalizer = np.max(iDOS.get_DOS()[1])
                print(f"DOING: width={WIDTH}", flush=True)
                normalizer = 1 #! delete
                TDOS = TwistedDOS(mesh_TDMs, len(GM_set), round(np.rad2deg(theta), 6), cutoff=cutoff, 
                                  width=WIDTH, kdim=args.dos, normalizer=normalizer, eigsys=eigsys)
                if eigsys is None:
                    eigsys = TDOS.get_eigsys()
                omegas, DOS = TDOS.get_DOS()
                TPLT = TwistedPlotter(round(np.rad2deg(theta), 6), omegas, DOS, mode_set, corner_kmags, cutoff=cutoff)
                TPLT.make_plot(outdir=outdir, name=name, filename=outname, width=WIDTH)
        print("Modes outputted.")

    succ("Successfully completed phonon mode analysis (Took %.3lfs)."%(time()-start))


