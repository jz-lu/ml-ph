from ____exit_with_error import exit_with_error
from __class_input import InputData
from __printers import print_start_msg, print_end_msg
from ___constants_misc import GENERAL_ERR_USAGE_MSG
from __begin_computation import begin_computation
from ___constants_config import BZE_OUT_FILENAME
from pymatgen.io.vasp.inputs import Poscar
import copy, sys, os
import numpy as np

# TODO items:
# Get the put directories of phonopy right with os.chdir() and adjust the movers in phonopy and the cleanup functions accordingly.

# Run only this file by hand. This is the only file that will hold "constants" outside of constants files, due to user inputting directory of input files.
# The purpose of this script is to parse the command line and pass it inot main to being the calculation pipeline.

start_time = print_start_msg()

# Import and check the command line arguments.
# First flag indicates the type of operation we are doing (0 for standard, int > 0 for configuration)
# Second flag is just the directory where the input files are stored, we will make it the root of our calculations
# Third flag is whether to do van der Waals forces (T/F) bool
# The remaining flags need to say which calculations to do.
cmdargs = list(copy.deepcopy(sys.argv))
if '-f' in cmdargs:
    filename = cmdargs[cmdargs.index('-f') + 1]
    print("Reading program inputs settings from file named '%s'.\n"%(filename))
    try:
        with open(filename) as f:
            cmdargs = tuple(f.read().splitlines())
    except:
        exit_with_error(GENERAL_ERR_USAGE_MSG + "\n\n\t" + filename + "is not a valid file.\n\n")
else:
    cmdargs = tuple(cmdargs[1:]) # Get rid of the name of the program in the first arg
    print("Reading program inputs settings from command line.\n")

user_input_settings = InputData(cmdargs)

bze_points = begin_computation(user_input_settings)
print("All computations complete! Returning any debug messages and exiting.")
if bze_points: # (b, z, e)
    print("\n\n[Debug] BZE POINTS:")
    for i in bze_points:
        print(i)
    print("\n")

# Exit success
print_end_msg(start_time)

