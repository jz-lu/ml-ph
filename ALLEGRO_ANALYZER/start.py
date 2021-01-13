from ____exit_with_error import exit_with_error
from __class_input import InputData
from __printers import print_start_msg, print_end_msg
from ___constants_misc import GENERAL_ERR_USAGE_MSG
from __begin_computation import begin_computation
import copy, sys, os

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
cmdargs = tuple(copy.deepcopy(sys.argv))
if '-f' in cmdargs:
    ind = cmdargs.index('-f')
    try:
        with open(cmdargs[ind]) as f:
            cmdargs = tuple(f.read().splitlines())
    except:
        exit_with_error(GENERAL_ERR_USAGE_MSG)

user_input_settings = InputData(cmdargs)

# Move home directory to the selected one.
os.chdir(user_input_settings.get_base_root_dir())

bze_points = begin_computation(user_input_settings)
if bze_points == None:
    print_end_msg(start_time)
    sys.exit(GENERAL_ERR_USAGE_MSG)

# Analyze results of the bze calculations
print_end_msg(start_time)

# TODO: result functions to plot e vs z, e vs b, z vs b. 
# Recall bze_points is a list of tuples (b, z, e) which is easy to unzip and plot using Zoe's plotting method.


