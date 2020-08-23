from ____exit_with_error import exit_with_error
from ___constants_config import *
from __class_input import InputData

from __printers import print_start_msg
from compute_properties import relax_solid
import copy
import sys
import multiprocessing as mp

# Run only this file by hand. This is the only file that will hold "constants" outside of constants files, due to user inputting directory of input files.
# The purpose of this script is to parse the command line and pass it inot main to being the calculation pipeline.

print_start_msg()

# Import and check the command line arguments. 
# First flag indicates the type of operation we are doing (0 for standard, int > 0 for configuration)
# Second flag is just the directory where the input files are stored, we will make it the root of our calculations
# Third flag is whether to do van der Waals forces (T/F) bool
# The remaining flags need to say which calculations to do.
cmdargs = tuple(copy.deepcopy(sys.argv))
user_input_settings = InputData(cmdargs)

# Build configuration sampling from user input
def config(user_input_settings):
    if user_input_settings.get_type_flag() == 0:
        print('Running standard computation.')
        relax_solid(user_input_settings)
    else:
        print('Running parallel computations over grid sample of configuration space, defaulted to layer 1 (z = 0)...')
        # Sample the grid here! Then pool them over to relax and run an iterator to get energy pairs for graphing
        # Num processes (in pool arg below) = number of grid points, i.e. (a, b, c) |-> a * b * c

        pool = mp.Pool()
