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

USAGE_ERR_MSG = 'Usage: python3 <DIR>/config_analyze.py -n <NUM SHIFTS> -d <I/O DIR FROM MAIN PROGRAM> -e <MIN ENERGY (eV)>'

# The DataOutput class expects a COB matrix and 
# a list of (config vector b, z-spacing, energy in eV), as well as a minimum energy in eV to shift by.
# This was all outputted in various files in the calculations, which need to be parsed.

# Parse cmdline args
cmdargs = list(copy.deepcopy(sys.argv))[1:]; i = 0; n = len(cmdargs)
BASE_ROOT = None; abs_min_energy = None; nshifts = None
while i < n:
    if cmdargs[i] == '-n':
        i += 1; nshifts = int(cmdargs[i]); i += 1
    elif cmdargs[i] == '-d':
        i += 1; BASE_ROOT = checkPath(cmdargs[i]); i += 1
        assert isdir(BASE_ROOT) # specify full path if not working properly
    elif cmdargs[i] == '-e':
        i += 1; abs_min_energy = float(cmdargs[i]); i += 1
    else:
        print(USAGE_ERR_MSG)
        sys.exit()
assert BASE_ROOT and abs_min_energy

print("== Configuration Analyzer Starting =="); start_time = time()
print("WD: %s, number of shifts: %d, minimum energy (eV): %.6lf."%(BASE_ROOT, nshifts, abs_min_energy))

# Collect COB matrix
data_dir = BASE_ROOT + checkPath(CONFIG_DATA_DIR)
print("Retrieving COB matrix from '%s'..."%data_dir)
cob = np.load(data_dir + COB_NPY_NAME)
print("COB matrix retrieved.")

# Collect shifts (including the 0 in the z direction)
b = [0]*nshifts
for i in range(nshifts):
    with open(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + SHIFT_NAME, 'r') as f:
        b[i] = tuple(map(float, f.read().splitlines()))
# Combine into (b, z, e) points and pass to DataOutput
bze = []

# Collect z-spacings
z = [0]*nshifts
for i in range(nshifts):
    relax_dir = BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) + checkPath(RELAXATION_DIR_NAME)
    z[i] = CarCollector.get_interlayer_spacing(relax_dir)

# Collect energies
print("Retrieving energies...")
e = [0]*nshifts
for i in range(nshifts):
    with open(BASE_ROOT + checkPath(CONFIG_SUBDIR_NAME + str(i)) 
                + checkPath(ANALYSIS_DIR_NAME) + checkPath(TOTAL_ENER_DIR_NAME) + TOT_ENERGIES_NAME) as f:
        e[i] = float(f.readline.split(' ')[-1])
print("Energies retrieved.")

bze = np.array(bze)
np.save(data_dir + 'bze', bze)
do = None
print("Parsing successful, passing to analyzer...")
if abs_min_energy is None:
    do = DataOutput(data_dir, bze, cob)
else:
    do = DataOutput(data_dir, bze, cob, abs_min_energy)
do.output_all_analysis()
print("Analyzer has finished running.")

print("== Configuration Analyzer Complete (Time: %.3lf s) =="%(time()-start_time))
