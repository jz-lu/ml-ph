# Main file to be executed by user
import sys
from ___constants_names import *
from ___constants_misc import *

# TODO: get a function for total energy and fermi energy
# TODO (s): future, in command line arguments add a setting of alreadyRelaxed so that if so we skip the relaxation and go to calculation
# TODO (s): get the automatic k-point line generation

# First check the command line arguments. They need to say which calculations to do.
calculation_list = sys.argv
if len(calculation_list) == 1:
        sys.exit(GENERAL_ERR_USAGE_MSG) # Need at least one calculation

calculation_list = list(dict.fromkeys(calculation_list)) # converts to dict, which removes duplicates, and back to list
del calculation_list[0]
for i in calculation_list:
    if i not in CMD_LINE_ARG_LIST:
        print(GENERAL_ERR_USAGE_MSG)
        sys.exit(BAD_INPUT_ERR_MSG)
    # Makes sure that we do not have invalid command line arguments before we begin postprocessing

# Check for the existence of a line kpoints file
# NOTE: remove this later if we can get the auto kpoints generated.

for i in calculation_list:
    if i == ENERGIES:
        # Generate the outcar and print total energies to a file here.




