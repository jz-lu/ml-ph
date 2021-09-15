"""
Script to analyze output of ml-ph 
VASP calculation for optimal interlayer spacing.
"""
from ___constants_names import (
    ANALYSIS_DIR_NAME, RELAXATION_DIR_NAME, 
    TOTAL_ENER_DIR_NAME, TOT_ENERGIES_NAME, 
    CONFIG_SUBDIR_NAME, CONFIG_DATA_DIR, 
    LIDXS_NPY_NAME
)
from __directory_searchers import findDirsinDir, checkPath
from __class_CarCollector import CarCollector
import matplotlib.pyplot as plt
import os, argparse
import numpy as np

def optim_interlayer_spacing(ROOT, plot=True):
    ROOT = checkPath(os.path.abspath(ROOT))
    data_dir = ROOT + checkPath(CONFIG_DATA_DIR)
    lidxs = np.load(data_dir + LIDXS_NPY_NAME)
    nshifts = len(findDirsinDir(ROOT, CONFIG_SUBDIR_NAME, searchType='start'))
    print(f"Sample size: {nshifts}")
    print("Retrieving z-spacings...")
    zspaces = np.zeros(nshifts)
    for i in range(nshifts):
        relax_dir = ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + checkPath(RELAXATION_DIR_NAME)
        zspaces[i] = CarCollector.get_interlayer_spacing(lidxs, DIR_RELAXATION=relax_dir)
    print("z-spacings retrieved.")

    # Collect energies
    print("Retrieving energies...")
    energies = np.zeros(nshifts)
    for i in range(nshifts):
        with open(ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) 
                    + checkPath(ANALYSIS_DIR_NAME) + checkPath(TOTAL_ENER_DIR_NAME) + TOT_ENERGIES_NAME) as f:
            energies[i] = float(f.readline().split(' ')[-1])
    print(f"Energies retrieved.")
    e = min(energies)
    idx = np.where(energies == e)[0]
    z = zspaces[idx]
    print("Cfg idx | interlayer spacings | energy:")
    print(np.vstack((np.arange(nshifts)+1, zspaces, (energies - e)*1000)))
    
    print(f"Index with minimum energy: {idx}, spacing (direct): {z}, energy: {'%.3lf'%e} eV")
    if plot:
        plt.clf(); fig, ax = plt.subplots()
        plt.title("Energy vs. interlayer spacing")
        plt.xlabel("Interlayer spacing (direct)"); plt.ylabel("Energy (eV)")
        ax.scatter(zspaces, energies, c='black')
        ax.plot(zspaces, energies, c='royalblue')
        fig.savefig(ROOT + "einter.png")
        plt.close(fig)
    return idx, z, e


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Optimize interlayer spacing for bilayer AA stacking")
    parser.add_argument("-d", "--dir", type=str, help="directory containing calculations", default='.')
    parser.add_argument("-n", "--noplot", action="store_true", help="stop plotting")
    args = parser.parse_args()
    
    assert os.path.isdir(args.dir), f"Directory {args.dir} does not exist"
    optim_interlayer_spacing(args.dir, plot=not args.noplot)
    

