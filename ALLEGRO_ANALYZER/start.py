import sys, os
from datetime import datetime
from pymatgen.io.vasp.inputs import Kpoints, Poscar, VaspInput
from pymatgen.io.vasp.outputs import Outcar

from ____exit_with_error import exit_with_error

from ___constants_misc import *
from ___constants_names import *
from ___constants_vasp import VASP_OUTFILE_LEN_THRESHOLD

from __directory_searchers import checkPath
from __dirModifications import mkdir
from __build_inputs import buildInitialInputs
from __run_vasp import run_vasp
from __get_energy_analysis import get_energy_analysis
from __postprocess_relaxation import postProcess_relaxation

print(PRGM_COOL_NAME)
print('\n')
print('Version:', PRGM_VERSION)
print(datetime.now(), '\n')

# Run only this file by hand. This is the only file that will hold "constants" outside of constants files, due to user inputting directory of input files.
# The purpose of this script is to parse the command line and pass it inot main to being the calculation pipeline.

needEnergies = False # Total energy requires no further processing after relaxation, so we treat it separately

# First check the command line arguments. 
# The first flag is just the directory where the input files are stored, we will make it the root of our calculations
# Second flag is whether to do van der Waals forces (T/F) bool
# The remaining flags need to say which calculations to do.
cmdargs = sys.argv

print('Command line arguments imported, now parsing... Arguments: ', cmdargs)

if len(cmdargs) <= 3:
        exit_with_error(GENERAL_ERR_USAGE_MSG) # Need at least one calculation

# In order of the flags, we store and check validity of the input arguments
del cmdargs[0] # First element is program name
# Get the root directory from the command line and check validity
ROOT = checkPath(cmdargs.pop(0))
if not os.path.isdir(ROOT):
    exit_with_error(ERR_BAD_INPUT_DIR)

# Check van der Waals flag following root dir specification above
do_vdW = cmdargs.pop(0)
if do_vdW == 'T' or do_vdW == 'True':
    do_vdW = True
    print('Turning on van der Waals settings...')
elif do_vdW == 'F' or do_vdW == 'False':
    do_vdW = False
else:
    print(ERR_INVALID_VDW_FLAG)
    exit_with_error(GENERAL_ERR_USAGE_MSG)

kpoints_is_gamma_centered = False
if cmdargs.pop(0).lower()[0] == 'g': # i.e. if user entered Gamma-centered kpoints
    kpoints_is_gamma_centered = True
    print('Adopting Gamma-centered sampling of the Brillouin zone.')
else:
    print('Adoping Monkhorst-Pack scheme of sampling the Brillouin zone.')

# At this point the arguments not popped are all calculation flags so we store it as calculation_list
calculation_list = cmdargs = list(dict.fromkeys(cmdargs)) # converts to dict, which removes duplicates, and back to list 
print('Arguments parsed. List of calculations to perform:', calculation_list)

enerIndex = None
for i in range(0, len(calculation_list)):
        if calculation_list[i] not in CMD_LINE_ARG_LIST:
            print(GENERAL_ERR_USAGE_MSG)
            exit_with_error(BAD_INPUT_ERR_MSG)
        # Makes sure that we do not have invalid command line arguments before we begin postprocessig
        if calculation_list[i] == ENERGIES:
            enerIndex = i
            needEnergies = True
if enerIndex != None:
    del calculation_list[i] # We will parse the total energy separately and before everything else since we can stop there if nothing else is needed.
# calculation_list is now fully parsed and ready to send to postprocessing

# Try to get the poscar, if fail then exit
print('Importing POSCAR from specified file directory %s...'%(ROOT))
try:
    poscar = Poscar.from_file(ROOT + POSCAR_UNIT_RELAXATION_NAME)
except Exception as err:
    print('Error:', err)
    exit_with_error(ERR_NO_POSCAR)

# Before we build the relaxation calculation inputs below, we need a sanity check as to whether the calculations can be done
# The only thing we need to check besides poscar validity (code will handle other input files if none are given
# is the existence of a line kpoints file if we want an electronic band structure calculation
kpoints_line = None
if (ELEBAND in calculation_list) or (PHBAND in calculation_list):
    try:
        print('Parsing imported KPOINTS line file for band structure calculations...')
        kpoints_line = Kpoints.from_file(ROOT + KPOINTS_LINE_NAME)
    except Exception as err:
        exit_with_error('Error in importing line KPOINTS file for requested electronic and/or phononic band structure calculations: ' + err)
# If we need the band calculations then get the line kpoints first (i.e. check to make sure it's valid, then might as well get it now)

# Is there a POTCAR given already? or should we generate a new one? Let's find out here.
potcarExists = os.path.isfile(ROOT + POTCAR_NAME)

# Build the necessary input files, returns vasp object
init_vasp_obj = buildInitialInputs(ROOT, do_vdW, poscar, kpoints_is_gamma_centered, potcarExists)

# We want to run the relaxation on a separate subfolder, so everything is organized
mkdir(RELAXATION_DIR_NAME, ROOT)
DIR_RELAXATION = checkPath(ROOT + RELAXATION_DIR_NAME)
print('Created new folder %s to store relaxation calculations.'%(DIR_RELAXATION))

# Similarly we want to run the analyses post-relaxation in a separate subfolder
mkdir(ANALYSIS_DIR_NAME, ROOT)
DIR_ANALYSIS = checkPath(ROOT + ANALYSIS_DIR_NAME)
print('Created new folder %s to store relaxation calculations.'%(DIR_ANALYSIS))

# Call the relaxation
print('Running VASP relaxation calculations...')
print(init_vasp_obj)
print(DIR_RELAXATION)
run_vasp(init_vasp_obj, DIR_RELAXATION)
print('VASP relaxation calculations complete.')

# Begin post-processing the relaxation calculations, handled by the postprocessing module
# Step 1: do the total energy
if needEnergies:
    outcar = Outcar(DIR_RELAXATION + OUTCAR_NAME)
    mkdir(TOTAL_ENER_DIR_NAME, DIR_ANALYSIS)
    energyPair = get_energy_analysis(checkPath(DIR_ANALYSIS + TOTAL_ENER_DIR_NAME), outcar) # energyPair = (total energy, fermi energy) tuple pair
if len(calculation_list) == 0:
    exit_with_error('Successfully completed total energy calculation, which was the only specified calculation.')

# Step 2, for ele/ph calculations, we transfer to postprocessing module
postProcess_relaxation(DIR_ANALYSIS, DIR_RELAXATION, init_vasp_obj, calculation_list, kpoints_line)