# Parses data output from main program and passes to analyzer
from ___constants_names import (
    CONFIG_SUBDIR_NAME, 
    CONFIG_DATA_DIR, 
    SHIFT_NAME, 
    COB_NPY_NAME, 
    RELAXATION_DIR_NAME, 
    ANALYSIS_DIR_NAME, 
    TOTAL_ENER_DIR_NAME, 
    TOT_ENERGIES_NAME
)
from __class_DataOutput import DataOutput
from __directory_searchers import checkPath
from __class_CarCollector import CarCollector
import numpy as np
import sys, copy
from os.path import isdir
from time import time
from ___helpers_parsing import succ, warn, err

USAGE_ERR_MSG = 'Usage: python3 <DIR>/config_analyze.py -n <NUM SHIFTS> -d <I/O DIR FROM MAIN PROGRAM> -e <MIN ENERGY (eV)>'

# The DataOutput class expects a COB matrix and 
# a list of (config vector b, z-spacing, energy in eV), as well as a minimum energy in eV to shift by.
# This was all outputted in various files in the calculations, which need to be parsed.

# Parse cmdline args
cmdargs = list(copy.deepcopy(sys.argv))[1:]; i = 0; n = len(cmdargs)
BASE_ROOT = None; abs_min_energy = None; nshifts = None
nlevel = 301
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
    else:
        warn(f"Unrecognized token '{cmdargs[i]}'")
        err(USAGE_ERR_MSG)
assert BASE_ROOT and nlevel > 0

print("== Configuration Analyzer Starting =="); start_time = time()
print("WD: %s, number of shifts: %d."%(BASE_ROOT, nshifts))
if not abs_min_energy:
    print("Using automatic minimum energy shift.")
else:
    print(f"Using minimum energy shift {abs_min_energy} eV.")

# Collect COB matrix
data_dir = BASE_ROOT + checkPath(CONFIG_DATA_DIR)
print("Retrieving COB matrix from '%s'..."%data_dir)
cob = np.load(data_dir + COB_NPY_NAME)
print("COB matrix retrieved.")

# Collect shifts (including the 0 in the z direction)
bshifts = [0]*nshifts
print("Retrieving shift coordinates...")
for i in range(nshifts):
    with open(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + SHIFT_NAME, 'r') as f:
        bshifts[i] = np.array(list(map(float, f.read().splitlines())))
print("Shift coordinates retrieved.")

# Collect z-spacings
print("Retrieving z-spacings...")
zspaces = [0]*nshifts
for i in range(nshifts):
    relax_dir = BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + checkPath(RELAXATION_DIR_NAME)
    zspaces[i] = CarCollector.get_interlayer_spacing(relax_dir)
print("z-spacings retrieved.")

# Collect energies
print("Retrieving energies...")
energies = [0]*nshifts
for i in range(nshifts):
    with open(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) 
                + checkPath(ANALYSIS_DIR_NAME) + checkPath(TOTAL_ENER_DIR_NAME) + TOT_ENERGIES_NAME) as f:
        energies[i] = float(f.readline().split(' ')[-1])
print("Energies retrieved.")

# Combine into (b, z, e) points and pass to DataOutput
bze = []
print("Combining data into bze-points data structure...")
for b, z, e in zip(bshifts, zspaces, energies):
    bze.append(np.array([b, z, e]))
print("Successfully combined.")

bze = np.array(bze)
np.save(data_dir + 'bze', bze)
do = None
print("Parsing successful, passing to analyzer...")
do = DataOutput(data_dir, bze, cob, abs_min_energy=abs_min_energy)
do.output_all_analysis(levels=nlevel)
print("Analyzer has finished running.")

succ("== Configuration Analyzer Complete (Took %.3lfs) =="%(time()-start_time))
