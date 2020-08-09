# Main file to be executed by user
import os
from pymatgen.io.vasp.inputs import VaspInput

from ____exit_with_error import exit_with_error

from ___constants_names import BATCH_FILE_PATH, VASP_RUN_RELAX_OUT_NAME, VASP_RUN_RELAX_ERR_NAME
from ___constants_misc import ERR_VASP_RUN_RELAX
from ___constants_vasp import VASP_OUTFILE_LEN_THRESHOLD

from __directory_searchers import checkPath

# TODO (s): future, in command line arguments add a setting of alreadyRelaxed so that if so we skip the relaxation and go to calculation
# TODO (s): get the automatic k-point line generation

def run_vasp(vaspObj, dirName):
    try:
        dirName = checkPath(dirName)
        vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=VASP_RUN_RELAX_OUT_NAME, err_file=VASP_RUN_RELAX_ERR_NAME)
        
        # Check the most obvious error in the run automatically
        if os.stat(dirName + VASP_RUN_RELAX_OUT_NAME).st_size < VASP_OUTFILE_LEN_THRESHOLD:
            exit_with_error('It is very likely that something went wrong in the VASP relaxation calculation as the .out file is unreasonably short. If that is not the case, then modify the source code threshold constant.')
        
        return
    except Exception as err:
        exit_with_error(ERR_VASP_RUN_RELAX + ' ' + err)



