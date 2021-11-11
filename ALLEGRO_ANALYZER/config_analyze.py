# Parses data output from main program and passes to analyzer
from ___constants_names import (
    CONFIG_SUBDIR_NAME, CONFIG_DATA_DIR, 
    SHIFT_NAME, 
    COB_NPY_NAME, LIDXS_NPY_NAME, 
    RELAXATION_DIR_NAME, ANALYSIS_DIR_NAME, 
    TOTAL_ENER_DIR_NAME, TOT_ENERGIES_NAME, 
    POSCAR_NAME, POSCAR_CONFIG_NAMEPRE
)
from ___constants_output import NPREDPTS, DEFAULT_CONTOUR_LEVELS
from __class_ConfigOutput import ConfigOutput, DSamplingOutput, FourierGSFE
from __directory_searchers import checkPath, findFilesInDir
from __class_CarCollector import CarCollector
from __class_Configuration import Configuration
from __class_PhonopyAPI import PhonopyAPI
from __class_PhononConfig import PhononConfig
import numpy as np; import numpy.linalg as LA
import sys, copy
from os.path import isdir
import os
from time import time
from math import sqrt
import itertools
import argparse
from pymatgen.io.vasp.inputs import Poscar
from ___helpers_parsing import succ, warn, err, update, greet

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Analysis of configuration space calculations")
    parser.add_argument("-d", "--dir", type=str, help="main directory", default='.')
    parser.add_argument("-e", "--ame", type=float, help="absolute minimum energy (default auto)")
    parser.add_argument("-l", "--nlevel", type=int, help="number of contour levels", default=DEFAULT_CONTOUR_LEVELS)
    parser.add_argument("--diag", action="store_true", help="sampling is on a diagonal cut")
    parser.add_argument("--nff", action="store_true", help="do not do GSFE Fourier fit")
    parser.add_argument("--fc", action="store_true", help="compute force constants")
    parser.add_argument("--int", action="store_true", help="force constants: interlayer only")
    parser.add_argument("--cfx", action="store_true", help="force constants: same Cartesian coordinate only")
    parser.add_argument("--fcf", nargs="+", help="force constants: fixed set of choices")
    parser.add_argument("--phcfg", type=int, help="phonon configuration calculation on shift (X)", default=0)
    parser.add_argument("--perr", action="store_true", help="plot percent error in GSFE prediction")
    parser.add_argument("nshifts", type=int, help='enter N^2 for NxN grid, or N for N diagonal cut', default=81)
    args = parser.parse_args()

    # Legacy compatibility
    BASE_ROOT = checkPath(os.path.abspath(args.dir))
    nshifts = args.nshifts
    abs_min_energy = args.ame
    nlevel = args.nlevel
    diag = args.diag
    ff = not args.nff
    fc = args.fc
    INTERLAYER_ONLY = args.int
    SAME_CART_ONLY = args.cfx
    fclist = np.array(list(map(int, args.fcf))) if args.fcf is not None else None
    phcfg = args.phcfg
    plot_perr = args.perr
    nplt = int(sqrt(len(fclist))) if fclist is not None else 4
    fc_fixed = (args.fcf is not None)

    # The ConfigOutput class expects a COB matrix and 
    # a list of (config vector b, z-spacing, energy in eV), as well as a minimum energy in eV to shift by.
    # This was all outputted in various files in the calculations, which need to be parsed.

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

    poscars_uc = sorted(findFilesInDir(BASE_ROOT, POSCAR_CONFIG_NAMEPRE, searchType='start')); npos = len(poscars_uc)
    name = '-'.join([Poscar.from_file(BASE_ROOT + pname).comment for pname in poscars_uc])
    print(f"Name: {name}")

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
    print(f"Energies retrieved.")

    # cob = Poscar.from_file('POSCAR').structure.lattice.matrix[:2,:2].T
    # import csv
    # with open('e.txt', 'r') as f:
    #     c = csv.reader(f)
    #     rows = [row for row in c][1:]
    #     bshifts = np.array([list(map(float, row[:2])) for row in rows])
    #     bshifts = bshifts @ LA.inv(cob).T
    #     zspaces = np.array([float(row[2]) for row in rows])
    #     energies = np.array([float(row[3]) for row in rows])
    #     name = 'blG'
    #     lidxs = []

    ff_pred = None; large_cfg = None; do = None; ph_list = None
    if fc or phcfg >= 0:
        print("Building list of phonopy objects...")
        ph_api = PhonopyAPI(BASE_ROOT, ctype='config')
        ncfg, ph_list = ph_api.nconfigs(), ph_api.inter_ph_list()
        assert ncfg == nshifts, f"Number of configurations found {ncfg} inconsistent with number entered {nshifts}"
    if phcfg >= 0:
        print("Analyzing phonons in configuration space...")
        phonon_config = PhononConfig(bshifts, cob, ph_list, data_dir)
        phcfg_poscar = Poscar.from_file(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + POSCAR_NAME)
        phonon_config.plot_mode_quiver(phcfg_poscar, shift=phcfg)
        phonon_config.plot_mode_config()
        print("Phonons in configuration space analysis completed.")
    if ff:
        print("Running 6-term Fourier fitting on GSFE")
        stype = f'{nshifts if diag else (int(sqrt(nshifts)), int(sqrt(nshifts)))}-{"diagonal" if diag else "grid"}'
        ff_gsfe = FourierGSFE(energies, bshifts, cob, sampling_type=stype)
        ff_gsfe.output_all_analysis(data_dir)
        if diag:
            large_cfg = Configuration.sample_line(npts=NPREDPTS, basis=Configuration.diagonal_basis_from_cob(cob), dump=False, as_mat=True)
        else:
            large_cfg = np.array(Configuration.sample_grid(grid=(81, 81, 1)))
        ff_pred = ff_gsfe.predict(large_cfg)
    if diag:
        pts = [0, nshifts//3, 2*nshifts//3]
        update(f"Parsing successful (special points: {pts}), passing to analyzer...")
        if fc:
            print("Building list of phonopy objects...")
            ph_api = PhonopyAPI(BASE_ROOT, ctype='config')
            ncfg, ph_list = ph_api.nconfigs(), ph_api.inter_ph_list()
            assert ncfg == nshifts, f"Number of configurations found {ncfg} inconsistent with number entered {nshifts}"
        do = DSamplingOutput(data_dir, nshifts, name, special_pts=pts, energies=energies, spacings=zspaces, ph_list=ph_list)
        do.output_all_analysis()
        if ff_pred is not None:
            pts = [0, NPREDPTS//3, 2*NPREDPTS//3]
            addendum = 1000*(np.array(energies)-min(energies))
            addendum = (np.linspace(0, 1, num=nshifts+1), np.append(addendum, addendum[0]))
            do_pred = DSamplingOutput(data_dir, NPREDPTS, name, special_pts=pts, energies=ff_pred, scaled=True, dump=False)
            do_pred.plot_energies(pfx='pred', tsfx='fitted', interp=False, scat=False, line=True, addendum=addendum)
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
        do = ConfigOutput(data_dir, bze, cob, name, ph_list=ph_list, abs_min_energy=abs_min_energy)
        do.output_all_analysis(levels=nlevel)
        if ff_pred is not None:
            print("Building fitted GSFE plot...")
            do.plot_e_vs_b(pfx='pred', tpfx='Fitted', energies=ff_pred, b=large_cfg)
            do.plot_diag_cut(ff_pred, large_cfg, pfx='cut')
            if plot_perr:
                bmat = input("Enter path for b matrix .npy file: ")
                while not os.path.isfile(bmat):
                    bmat = input("Invalid path, try again: ")
                bmat = np.load(bmat)
                earr = input("Enter path for actual energies .npy file: ")
                while not os.path.isfile(earr):
                    earr = input("Invalid path, try again: ")
                earr = np.load(earr)
                ff_gsfe.plot_percent_error(bmat, earr)
                print("Successfully plotted percent error")
    if fc: # plot force constants
        update("Calculating force constants...")
        print("Getting POSCAR...")
        p = Poscar.from_file(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(0)) + POSCAR_NAME)
        n_at = len(p.structure.species); npairs = n_at * (n_at-1) // 2; at_idxs = np.arange(n_at)
        assert n_at == len(lidxs), f"Inconsistent number of atoms {n_at} with associated layer indices {lidxs}"
        assert nplt > 0, f"Invalid number of force constants to plot {nplt}"
        atomic_pairs = None; lidx_pairs = None; pfx = ''
        if fc_fixed:
            print("Using fixed pairs...")
            atomic_pairs = np.array([fclist[:2]])
            lidx_pairs = np.array([fclist[2:4]])
            cart_pairs = np.array([fclist[4:]])
            do.plot_forces(atomic_pairs, lidx_pairs, cart_pairs, p, pfx=f'fixed{"".join(list(map(str, fclist)))}')
        else:
            print("Sampling pairs...")
            assert nplt <= npairs, f"Number of force constants to plot {nplt} is too large, max={npairs}"
            atomic_pairs = np.array(list(itertools.combinations(at_idxs, 2)))
            lidx_pairs = np.array(list(itertools.combinations(lidxs, 2)))
            choice_idxs = np.random.choice(np.arange(len(atomic_pairs)), size=nplt, replace=False)
            print(f"Using pairs with indices: {choice_idxs}")
            atomic_pairs = atomic_pairs[choice_idxs]; lidx_pairs = lidx_pairs[choice_idxs]
            cart_pairs = np.stack((np.random.choice([0,1,2], nplt), np.random.choice([0,1,2], nplt)), axis=1)
            print("Plotting forces (no constraints)...")
            do.plot_forces(atomic_pairs, lidx_pairs, cart_pairs, p, pfx=pfx)
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
                cart_pairs = np.stack((np.random.choice([0,1,2], nplt), np.random.choice([0,1,2], nplt)), axis=1)
                print("Plotting forces (constraints: interlayer only)...")
                do.plot_forces(atomic_pairs, lidx_pairs, cart_pairs, p, pfx=pfx)
            if SAME_CART_ONLY:
                pfx += ('_' if INTERLAYER_ONLY else '') + 'cfix'
                atomic_pairs = np.array([[i, j] for i in at_idxs[lidxs==1] for j in at_idxs[lidxs==2]])
                np.random.shuffle(atomic_pairs)
                atomic_pairs = atomic_pairs[:nplt]
                carts = np.random.choice([0,1,2], nplt)
                cart_pairs = np.stack((carts, carts), axis=1)
                print("Plotting forces (constraints: interlayer only, same Cartesian coordinate)...")
                do.plot_forces(atomic_pairs, lidx_pairs, cart_pairs, p, pfx=pfx)
            
    print("Analyzer has finished running.")
    succ("== Configuration Analyzer Complete (Took %.3lfs) =="%(time()-start_time))

else:
    err(f"Error: '{__name__}' has nothing to import")