"""
Script to analyze the output of ml-ph 
VASP calculation for optimal lattice constant.
"""
from ___constants_names import (
    ANALYSIS_DIR_NAME, RELAXATION_DIR_NAME, 
    TOTAL_ENER_DIR_NAME, TOT_ENERGIES_NAME, 
    CONFIG_SUBDIR_NAME, CONFIG_DATA_DIR, 
    POSCAR_NAME
)
from __directory_searchers import findDirsinDir, checkPath
from pymatgen.io.vasp.inputs import Poscar
import matplotlib.pyplot as plt
import os, argparse
import numpy as np
import numpy.linalg as LA

def optim_lc(ROOT, plot=True):
    ROOT = checkPath(os.path.abspath(ROOT))
    data_dir = ROOT + checkPath(CONFIG_DATA_DIR)
    nshifts = len(findDirsinDir(ROOT, CONFIG_SUBDIR_NAME, searchType='start'))
    print(f"Sample size: {nshifts}")
    print("Retrieving lattice constants...")
    lcs = [0]*nshifts
    for i in range(nshifts):
        relax_dir = ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + checkPath(RELAXATION_DIR_NAME)
        lcs[i] = LA.norm(Poscar.from_file(relax_dir + POSCAR_NAME).structure.lattice.matrix[0])
    print("Lattice constants retrieved.")

    # Collect energies
    print("Retrieving energies...")
    energies = [0]*nshifts
    for i in range(nshifts):
        with open(ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) 
                    + checkPath(ANALYSIS_DIR_NAME) + checkPath(TOTAL_ENER_DIR_NAME) + TOT_ENERGIES_NAME) as f:
            energies[i] = float(f.readline().split(' ')[-1])
    print(f"Energies retrieved.")
    e = min(energies)
    idx = energies.index(min(energies))
    lc_star = lcs[idx]
    print(f"Index with minimum energy: {idx}, spacing (direct): {lc_star}, energy: {'%.3lf'%e} eV")
    if plot:
        plt.clf(); fig, ax = plt.subplots()
        plt.title("Energy vs. lattice constant")
        plt.xlabel(r"Lattice constants ($\AA$)"); plt.ylabel("Energy (eV)")
        ax.scatter(lcs, energies, c='black')
        ax.plot(lcs, energies, c='royalblue')
        fig.savefig(ROOT + "e_vs_lc.png")
        plt.close(fig)
    return idx, lc_star, e


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Optimize lattice constant for material calculation")
    parser.add_argument("-d", "--dir", type=str, help="directory containing calculations", default='.')
    parser.add_argument("-n", "--noplot", action="store_true", help="stop plotting")
    args = parser.parse_args()
    
    assert os.path.isdir(args.dir), f"Directory {args.dir} does not exist"
    optim_lc(args.dir, plot=not args.noplot)
    

