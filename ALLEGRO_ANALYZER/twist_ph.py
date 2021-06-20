# Compute the phonon modes for a twisted cell
import numpy as np
from phonopy import Phonopy
import phonopy
from __class_DM import TwistedDM, InterlayerDM, MonolayerDM
from __class_PhonopyAPI import PhonopyAPI
from bzsampler import BZSampler
from ___constants_names import (
    SPOSCAR_NAME, PH_FORCE_SETS_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME
)
from __directory_searchers import checkPath, findDirsinDir
from __dirModifications import build_dir
from ___helpers_parsing import greet, succ, warn, err, is_flag, check_not_flag
import os, sys

# Sample moire G vectors (GM) and k-points along IBZ


if __name__ == '__main__':
    USAGE_MSG = f"Usage: python3 {sys.argv[0]} -dir <directory given in twist calculation>"
    greet("Building twisted crystal phonon modes...")
    warn("Note: you must have already ran twisted calculation in given directory and retrieved forces from phonopy")
    args = sys.argv[1:]; i = 0; n = len(args)
    indir = None
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-dir':
            i += 1; check_not_flag(args[i])
            indir = args[i]
            if args[i] in ['.', './', '..', '../'] or args[i][0] == '.':
                warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
            i += 1
        else:
            warn(f'Warning: unknown flag "{args[i]} ignored')
            i += 1
    if not indir:
        err(USAGE_MSG)
    
    ph_api = PhonopyAPI(indir) # retreive phonopy objects
    




