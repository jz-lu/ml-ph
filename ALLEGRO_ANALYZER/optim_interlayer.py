from ___constants_names import (
    ANALYSIS_DIR_NAME, 
    TOTAL_ENER_DIR_NAME, TOT_ENERGIES_NAME, 
    CONFIG_SUBDIR_NAME
)
import matplotlib.pyplot as plt
import os, argparse

def optim_interlayer_spacing(ROOT, plot=True):
    return

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Optimize interlayer spacing for bilayer AA stacking")
    parser.add_argument("-d", "--dir", type=str, help="directory containing calculations", default='.')
    parser.add_argument("-n", "--noplot", action="store_true", help="stop plotting")
    args = parser.parse_args()
    
    assert os.path.isdir(args.dir), f"Directory {args.dir} does not exist"
    optim_interlayer_spacing(args.dir, plot=not args.noplot)
    

