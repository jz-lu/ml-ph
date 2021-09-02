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
    DEFAULT_PH_BAND_PLOT_NAME, 
    ANGLE_SAMPLE_INAME, 
    MODE_TNSR_ONAME, ANGLE_SAMPLE_ONAME, GAMMA_IDX_ONAME, K_MAGS_ONAME
)
from ___constants_compute import DEFAULT_NUM_LAYERS
from __directory_searchers import checkPath, findDirsinDir, findFilesInDir
from __dirModifications import build_dir
from __class_Configuration import Configuration
from __class_RelaxerAPI import RelaxerAPI
from __class_ForceInterp import ForceInterp, FourierForceInterp
from __class_PhononConfig import TwistedRealspacePhonon
from ___helpers_parsing import greet, update, succ, warn, err, is_flag, check_not_flag
import os, sys
from math import pi

RELAX_FOURIER_INTERP = 1
RELAX_SPLINE_INTERP = 2
    
if __name__ == '__main__':
    start = time()
    USAGE_MSG = f"Usage: python3 {sys.argv[0]} -deg <twist angle (deg)> -name <solid name> -cut <frequency cutoff> -dir <base I/O dir> -o <output dir>"
    args = sys.argv[1:]; i = 0; n = len(args)
    indir = '.'; outdir = '.'; theta = None; name = None; outname = DEFAULT_PH_BAND_PLOT_NAME; cutoff = None
    plot_intra = False; relax = False; realspace = False; force_sum = False; multirelax = False
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
        elif args[i] == '-name':
            i += 1; check_not_flag(args[i]); name = args[i]; i += 1
        elif args[i] == '-cut':
            i += 1; check_not_flag(args[i]); cutoff = float(args[i]); i += 1
        elif args[i] == '-fout':
            i += 1; check_not_flag(args[i]); outname = args[i]; i += 1
        elif args[i] == '--intra':
            i += 1; plot_intra = True
        elif args[i] == '--relax' or args[i] == '-r':
            i += 1
            relax = RELAX_SPLINE_INTERP if args[i] in ['r', 'real', 'R', 'REAL'] else RELAX_FOURIER_INTERP
            i += 1 
        elif args[i] == '--mr':
            i += 1
            multirelax = True
        elif args[i] == '--real':
            i += 1; realspace = True
        elif args[i] == '--fsum':
            i += 1; force_sum = True
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
        assert not relax, f"Cannot multirelax and relax"
        assert os.path.isfile(indir + ANGLE_SAMPLE_INAME), f"Must provide input file {ANGLE_SAMPLE_INAME}"
    print(f"WD: {indir}, Output to: {outdir}")
    outname = 'd' + str(round(np.rad2deg(theta))) + "_" + outname
    print("Building twisted crystal phonon modes...")
    warn("Note: you must have already ran twisted calculation in given directory and retrieved forces from phonopy")

    ldir = "/Users/jonathanlu/Documents/lukas/INTRA_FORCE_C/"
    greet("Working on intralayer components...")
    print("Searching for POSCAR inputs...")
    poscars_uc = [Poscar.from_file(ldir + "POSCAR_unit") for _ in range(2)]
    print("POSCAR inputs found.")

    super_dim = [5,5,1]
    per_layer_at_idxs = np.array([[0, 1, 2, 6, 7, 8, 9, 10, 11], [3, 4, 5, 12, 13, 14, 15, 16, 17]])
    ml_ph_list = []
    os.chdir(ldir)
    for _ in range(2):
            ml_ph_list.append(phonopy.load(supercell_matrix=[1, 1, 1],
                                            primitive_matrix='auto',
                                            unitcell_filename="./POSCAR",
                                            force_constants_filename="./FORCE_CONSTANTS"))
            ml_ph_list[-1]._log_level = 0
    os.chdir(outdir)

    update("Sampling G and k sets...")
    bzsamples = get_bz_sample(theta, poscars_uc[0], outdir, make_plot=True, super_dim=super_dim)
    _, GM_set = bzsamples.get_GM_set(); _, G0_set = bzsamples.get_G0_set()
    Gamma_idx = bzsamples.get_Gamma_idx()
    k_set, k_mags = bzsamples.get_kpts(); corner_kmags = bzsamples.get_corner_kmags()
    k0_set, k0_mags = bzsamples.get_kpts0(); corner_k0mags = bzsamples.get_corner_kmags0()
    np.save(outdir + "km.npy", k_mags)
    np.save(outdir + "k.npy", k_set)
    np.save(outdir + "ck.npy", corner_kmags)
    print("Sampling complete.")
    os.chdir(outdir)

    update("Constructing intralayer dynamical matrix objects...")
    MLDMs = [MonolayerDM(uc, None, ph, GM_set, \
                         G0_set, k_set, Gamma_idx, k_mags=k_mags) \
                         for uc, ph in zip(poscars_uc, ml_ph_list)]
    if plot_intra:
        print("Plotting one intralayer component...")
        print("Plotting in pristine space...")
        MLDMs[0].plot_pristine_band(k0_set, k0_mags, corner_k0mags, outdir=outdir)
        print("Plotting in moire space..."); MLDMs[0].plot_sampled_l0()
        MLDMs[0].plot_band(k_mags, corner_kmags, name=name, outdir=outdir, filename='intra_'+outname, cutoff=cutoff)
    if force_sum:
        for i, MLDM in enumerate(MLDMs):
            print(f"Computing layer {i} force sums...")
            MLDM.print_force_sum()
    print("Intralayer DM objects constructed.")

    greet("Working on interlayer components...")
    print("Importing shift vectors...")
    ldir = "/Users/jonathanlu/Documents/lukas/INTER_FORCE_C/"
    os.chdir(ldir)
    nshift = 100
    gridsz = int(sqrt(nshift)); assert gridsz**2 == nshift, f"Number of shifts {nshift} must be a perfect square"
    b_set = np.load("/Users/jonathanlu/Documents/lukas/b.npy")
    b_set = b_set @ poscars_uc[0].structure.lattice.matrix[:2,:2] # make Cartesian
    print("Shift vectors imported.")
        
    print("Getting interlayer phonopy objects from API...")
    config_ph_list = []
    for i in range(10):
        for j in range(10):
            config_ph_list.append(phonopy.load(supercell_matrix=[1, 1, 1],
                                    primitive_matrix='auto',
                                    unitcell_filename=f"./POSCAR_{i}_{j}",
                                    force_constants_filename=f"./FORCE_CONSTANTS_{i}_{j}"))
    print("Phonopy objects retrieved.")

    print("Note: Using GM sampling set from intralayer calculations.")
    print("Constructing interlayer dynamical matrix objects...")
    ILDM = None
    bl_M = np.array([95.94, 95.94, 32.065, 32.065, 32.065, 32.065])
    print(f"Bilayer masses: {bl_M}")
    ILDM =  InterlayerDM(per_layer_at_idxs, bl_M, 
                            b_set, 
                            bzsamples.get_kpts0()[0], 
                            GM_set, G0_set, 
                            [p.structure.species for p in poscars_uc], 
                            ph_list=config_ph_list, 
                            force_matrices=None)
    A0 = poscars_uc[0].structure.lattice.matrix[:2,:2].T
    G0_mat = 2*pi*LA.inv(A0).T
    ILDM.plot_pristine_band(G0_mat, k0_set, k0_mags, corner_k0mags, outdir=outdir)
    succ("yay")
    sys.exit()
    # evals = LA.eigvals(ILDM.get_DM())
    # signs = (-1)*(evals < 0) + (evals > 0) # pull negative sign out of square root to plot imaginary frequencies
    # modes_k = signs * np.sqrt(np.abs(evals)) * (15.633302*33.356) # eV/Angs^2 -> THz ~ 15.633302; THz -> cm^-1 ~ 33.356
    # print("Interlayer modes (k-independent):", modes_k[modes_k != 0])
    print("Interlayer DM objects constructed.")


    print("Combining into a single twisted dynamical matrix object...")
    TDM = TwistedDM(MLDMs[0], MLDMs[1], ILDM, k_mags, [p.structure.species for p in poscars_uc], Gamma_idx)

    # TDM.plot_band(corner_kmags, np.rad2deg(theta), outdir=outdir, name=name, filename="nosum"+outname, cutoff=cutoff)
    # TDM.modes_built = False
    TDM.Lukas_sum_rule()
    print("Twisted dynamical matrix object constructed.")

    print(f"Diagonalizing and outputting modes with corners {corner_kmags}...")
    TDM.plot_band(corner_kmags, np.rad2deg(theta), outdir=outdir, name=name, filename=outname, cutoff=cutoff)
    print("Modes outputted.")

    succ("Successfully completed phonon mode analysis (Took %.3lfs)."%(time()-start))


