# Parses data output from main program and passes to analyzer
from ___constants_names import (
    CONFIG_SUBDIR_NAME, CONFIG_DATA_DIR, 
    SHIFT_NAME, 
    COB_NPY_NAME, LIDXS_NPY_NAME, 
    RELAXATION_DIR_NAME, ANALYSIS_DIR_NAME, 
    TOTAL_ENER_DIR_NAME, TOT_ENERGIES_NAME, 
    POSCAR_NAME
)
from ___constants_output import NPREDPTS, DEFAULT_CONTOUR_LEVELS
from __class_ConfigOutput import ConfigOutput, DSamplingOutput, FourierGSFE
from __directory_searchers import checkPath
from __class_CarCollector import CarCollector
from __class_Configuration import Configuration
from __class_PhonopyAPI import PhonopyAPI
import numpy as np
import sys, copy
from os.path import isdir
import os
from time import time
from math import sqrt
import itertools
from pymatgen.io.vasp.inputs import Poscar
from ___helpers_parsing import succ, warn, err, update, greet
INTERLAYER_ONLY = False
SAME_CART_ONLY = False

if __name__ == '__main__':
    USAGE_ERR_MSG = 'Usage: python3 <DIR>/config_analyze.py -n <NUM SHIFTS> -d <I/O DIR FROM MAIN PROGRAM> -e <MIN ENERGY (eV)> (optional: --diag --nff --fc)'

    # The ConfigOutput class expects a COB matrix and 
    # a list of (config vector b, z-spacing, energy in eV), as well as a minimum energy in eV to shift by.
    # This was all outputted in various files in the calculations, which need to be parsed.

    # Parse cmdline args
    cmdargs = list(copy.deepcopy(sys.argv))[1:]; i = 0; n = len(cmdargs)
    BASE_ROOT = '.'; abs_min_energy = None; nshifts = None; diag = False; ff = True; fc = False; nplt = 4
    nlevel = DEFAULT_CONTOUR_LEVELS
    while i < n:
        if cmdargs[i] == '-n':
            i += 1; nshifts = int(cmdargs[i]); i += 1
        elif cmdargs[i] == '-d':
            i += 1; BASE_ROOT = checkPath(cmdargs[i]); i += 1
            assert isdir(BASE_ROOT) # specify full path if not working properly
        elif cmdargs[i] == '-e':
            i += 1; abs_min_energy = float(cmdargs[i]); i += 1
        elif cmdargs[i] == '-l':
            i += 1; nlevel = int(cmdargs[i]); i += 1
        elif cmdargs[i] == '--diag':
            diag = True; i += 1
        elif cmdargs[i] == '--nff':
            ff = False; i += 1
        elif cmdargs[i] == '--int':
            INTERLAYER_ONLY = True; i += 1
        elif cmdargs[i] == '--cfx':
            SAME_CART_ONLY  = True; i += 1
        elif cmdargs[i] == '--fc':
            i += 1; nplt = int(cmdargs[i]); i +=1; fc = True
        elif cmdargs[i] == '--usage':
            print(USAGE_ERR_MSG)
            sys.exit(0)
        else:
            warn(f"Unrecognized token '{cmdargs[i]}'")
            err(USAGE_ERR_MSG)
    assert BASE_ROOT and nlevel > 0 and nplt > 0
    BASE_ROOT = checkPath(os.path.abspath(BASE_ROOT))
    if int(sqrt(nplt))**2 != nplt:
        print(f"Rounded number of forces {nplt} down to nearest perfect square {int(sqrt(nplt))}")
        nplt = int(sqrt(nplt))

    greet("== Configuration Analyzer Starting =="); start_time = time()
    update("WD: %s, number of shifts: %d."%(BASE_ROOT, nshifts))
    if not abs_min_energy:
        print("Using automatic minimum energy shift.")
    else:
        print(f"Using minimum energy shift {abs_min_energy} eV.")
    print(f"Analysis type: {'along diagonal' if diag else 'grid configurations'}")

    # Collect COB matrix and layer indices
    data_dir = BASE_ROOT + checkPath(CONFIG_DATA_DIR)
    print("Retrieving COB matrix and layer indices from '%s'..."%data_dir)
    cob = np.load(data_dir + COB_NPY_NAME)
    print("COB matrix retrieved.")
    lidxs = np.load(data_dir + LIDXS_NPY_NAME)
    print("Layer indices retrieved.")

    # Collect shifts (including the 0 in the z direction)
    bshifts = [0]*nshifts
    print("Retrieving shift coordinates...")
    for i in range(nshifts):
        with open(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + SHIFT_NAME, 'r') as f:
            bshifts[i] = list(map(float, f.read().splitlines()))
    bshifts = np.array(bshifts)
    np.save(data_dir + 'bdirect.npy', bshifts)
    np.save(data_dir + 'bcart.npy', (cob @ bshifts[:,:2].T).T)
    print("Shift coordinates retrieved.")

    # Collect z-spacings
    print("Retrieving z-spacings...")
    zspaces = [0]*nshifts
    for i in range(nshifts):
        relax_dir = BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + checkPath(RELAXATION_DIR_NAME)
        zspaces[i] = CarCollector.get_interlayer_spacing(lidxs, DIR_RELAXATION=relax_dir)
    print("z-spacings retrieved.")

    # Collect energies
    print("Retrieving energies...")
    energies = [0]*nshifts
    for i in range(nshifts):
        with open(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) 
                    + checkPath(ANALYSIS_DIR_NAME) + checkPath(TOTAL_ENER_DIR_NAME) + TOT_ENERGIES_NAME) as f:
            energies[i] = float(f.readline().split(' ')[-1])
    print(f"Energies retrieved: {energies}")

    ff_pred = None; large_cfg = None; do = None; ph_list = None
    if fc:
        print("Building list of phonopy objects...")
        ph_api = PhonopyAPI(BASE_ROOT, ctype='config')
        ncfg, ph_list = ph_api.nconfigs(), ph_api.inter_ph_list()
        assert ncfg == nshifts, f"Number of configurations found {ncfg} inconsistent with number entered {nshifts}"
    if ff:
        print("Running 6-term Fourier fitting on GSFE")
        stype = f'{nshifts if diag else (int(sqrt(nshifts)), int(sqrt(nshifts)))}-{"diagonal" if diag else "grid"}'
        ff_gsfe = FourierGSFE(energies, bshifts, cob, sampling_type=stype)
        ff_gsfe.output_all_analysis(data_dir)
        if diag:
            large_cfg = Configuration.sample_line(npts=NPREDPTS, basis=Configuration.diagonal_basis_from_cob(cob), dump=False, as_mat=True)
        else:
            large_cfg = np.array(Configuration.sample_grid(grid=(101, 101, 1)))
        ff_pred = ff_gsfe.predict(large_cfg)
    if diag:
        pts = [0, nshifts//3, 2*nshifts//3]
        update(f"Parsing successful (special points: {pts}), passing to analyzer...")
        if fc:
            print("Building list of phonopy objects...")
            ph_api = PhonopyAPI(BASE_ROOT, ctype='config')
            ncfg, ph_list = ph_api.nconfigs(), ph_api.inter_ph_list()
            assert ncfg == nshifts, f"Number of configurations found {ncfg} inconsistent with number entered {nshifts}"
        do = DSamplingOutput(data_dir, nshifts, special_pts=pts, energies=energies, spacings=zspaces, ph_list=ph_list)
        do.output_all_analysis()
        if ff_pred is not None:
            pts = [0, NPREDPTS//3, 2*NPREDPTS//3]
            do_pred = DSamplingOutput(data_dir, NPREDPTS, special_pts=pts, energies=ff_pred, scaled=True, dump=False)
            do_pred.plot_energies(pfx='pred', tsfx='fitted', interp=False, scat=False, line=True)
    else:
        # Combine into (b, z, e) points and pass to ConfigOutput
        bze = []
        update("Combining data into bze-points data structure...")
        for b, z, e in zip(bshifts, zspaces, energies):
            bze.append(np.array([b, z, e]))
        print("Successfully combined.")
        bze = np.array(bze)
        np.save(data_dir + 'bze', bze)
        print("Parsing successful, passing to analyzer...")
        do = ConfigOutput(data_dir, bze, cob, ph_list=ph_list, abs_min_energy=abs_min_energy)
        do.output_all_analysis(levels=nlevel)
        if ff_pred is not None:
            print("Building fitted GSFE plot...")
            do.plot_e_vs_b(pfx='pred', tpfx='Fitted', energies=ff_pred, b=large_cfg)
    if fc: # plot force constants
        assert nplt > 0, f"Invalid number of force constants to plot {nplt}"
        update("Calculating force constants...")
        print("Getting POSCAR...")
        p = Poscar.from_file(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(0)) + POSCAR_NAME)
        n_at = len(p.structure.species); npairs = n_at * (n_at-1) // 2; at_idxs = np.arange(n_at)
        assert n_at == len(lidxs), f"Inconsistent number of atoms {n_at} with associated layer indices {lidxs}"
        print("Sampling pairs...")
        atomic_pairs = None; lidx_pairs = None; pfx = ''
        if INTERLAYER_ONLY:
            print("Using constraint: interlayer forces only")
            n_at //= 2; npairs = (n_at * (n_at-1) // 2)**2
            assert nplt <= npairs, f"Number of force constants to plot {nplt} is too large, max={npairs}"
            atomic_pairs = np.array([[i, j] for i in at_idxs[lidxs==1] for j in at_idxs[lidxs==2]])
            np.random.shuffle(atomic_pairs)
            atomic_pairs = atomic_pairs[:nplt]
            lidx_pairs = np.array([[1,2]]*nplt)
            pfx += 'inter'
            print(f"Atomic pairs ({atomic_pairs.shape}):\n{atomic_pairs}\nlidx_pairs ({lidx_pairs.shape}): \n{lidx_pairs}")
        else:
            assert nplt <= npairs, f"Number of force constants to plot {nplt} is too large, max={npairs}"
            atomic_pairs = np.array(list(itertools.combinations(at_idxs, 2)))
            lidx_pairs = np.array(list(itertools.combinations(lidxs, 2)))
            choice_idxs = np.random.choice(np.arange(len(atomic_pairs)), size=nplt, replace=False)
            print(f"Using pairs with indices: {choice_idxs}")
            atomic_pairs = atomic_pairs[choice_idxs]; lidx_pairs = lidx_pairs[choice_idxs]
        cart_pairs = np.stack((np.random.choice([0,1,2], nplt), np.random.choice([0,1,2], nplt)), axis=1)
        if SAME_CART_ONLY:
            pfx += ('_' if INTERLAYER_ONLY else '') + 'cfix'
            carts = np.random.choice([0,1,2], nplt)
            cart_pairs = np.stack((carts, carts), axis=1)
        print("Plotting forces...")
        do.plot_forces(atomic_pairs, lidx_pairs, cart_pairs, p, pfx=pfx)
            
    print("Analyzer has finished running.")
    succ("== Configuration Analyzer Complete (Took %.3lfs) =="%(time()-start_time))

else:
    err("Error: script has nothing to import")