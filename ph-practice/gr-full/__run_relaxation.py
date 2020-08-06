# Main file to be executed by user
import sys
import os
from pymatgen.io.vasp.inputs import VaspInput

from ___constants_names import BATCH_FILE_PATH, VASP_RUN_RELAX_OUT_NAME, VASP_RUN_RELAX_ERR_NAME
from ___constants_misc import ERR_VASP_RUN_RELAX

from __directory_searchers import checkPath

# TODO: get a function for total energy and fermi energy
# TODO (s): future, in command line arguments add a setting of alreadyRelaxed so that if so we skip the relaxation and go to calculation
# TODO (s): get the automatic k-point line generation

def run_relaxation(vaspObj, dirName):
    try:
        dirName = checkPath(dirName)
        vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=VASP_RUN_RELAX_OUT_NAME, err_file=VASP_RUN_RELAX_ERR_NAME)
        return dirName + VASP_RUN_RELAX_OUT_NAME # for catching huge run errors
    except Exception as err:
        sys.exit(ERR_VASP_RUN_RELAX, err)




