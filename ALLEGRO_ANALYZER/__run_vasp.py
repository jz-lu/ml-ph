# Main file to be executed by user
import os
from pymatgen.io.vasp.inputs import VaspInput, Poscar
from ____exit_with_error import exit_with_error
from ____debug import DEBUG_NOTICE_MSG, DEBUGGING
from ___constants_names import (
    BATCH_FILE_PATH, 
    VASP_RUN_OUT_NAME, 
    VASP_RUN_ERR_NAME, 
    CHGCAR_NAME, 
    CONTCAR_NAME, 
    POSCAR_NAME, 
    INCAR_RELAXATION_NAME,
    KPOINTS_NAME
)
from ___constants_misc import ERR_VASP_RUN_RELAX, ERR_VASP_NOT_CONVERGED
from ___constants_vasp import VASP_OUTFILE_LEN_THRESHOLD, VASP_MAX_CONVERGENCE_ATTEMPTS
from __class_CarCollector import CarCollector
from __directory_searchers import checkPath
from __check_convergence import check_if_not_converged

def run_vasp(vaspObj, dirName, predefined_chgcar=None, run_type='relax'):
    if DEBUGGING:
        print(DEBUG_NOTICE_MSG)
        return
    dirName = checkPath(dirName)

    # pymatgen doesn't handle other files well. We will manually print chgcar in the outdir so VASP will take it when it runs.
    if predefined_chgcar != None:
        predefined_chgcar.write_file(dirName + CHGCAR_NAME)
    
    outfile_name = VASP_RUN_OUT_NAME%(run_type)
    errfile_name = VASP_RUN_ERR_NAME%(run_type)

    # No relaxation -> no need to check for ionic convergence
    if run_type != 'relax':
        vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=outfile_name, err_file=errfile_name)

    # Relaxation -> need to make sure it converges
    else:
        loopCounter = 1 # Counts number of times we try to relax. We stop at VASP_MAX_CONVERGENCE_ATTEMPTS

        # Run vasp and check convergence
        vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=outfile_name, err_file=errfile_name)
        not_converged = check_if_not_converged(dirName, outfile_name)
        if not not_converged:
            print('Relaxation has converged.')

        while not_converged:
            print('Relaxation has not converged. Proceeding to rerun relaxation with updated ionic positions...')
            print('Loop counter: ' + str(loopCounter))
            loopCounter += 1
            # Loop counter and exit function
            if loopCounter > VASP_MAX_CONVERGENCE_ATTEMPTS:
                # exit_with_error(ERR_VASP_NOT_CONVERGED)
                print(ERR_VASP_NOT_CONVERGED)
                break
            
            new_poscar = Poscar.from_file(dirName + CONTCAR_NAME)
            vaspObj['POSCAR'] = new_poscar # Update POSCAR to be the new CONTCAR
            print('CONTCAR copied to POSCAR. Rerunning relaxation...')
            # No need to delete anything, the files will be overwritten. Just start the calculation.
            vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=outfile_name, err_file=errfile_name)
            print('Run complete.')
            not_converged = check_if_not_converged(dirName, outfile_name)
            if not not_converged:
                print('Relaxation has converged.')
    
    # Check the most obvious error in the run automatically
    if os.stat(dirName + outfile_name).st_size < VASP_OUTFILE_LEN_THRESHOLD:
        exit_with_error('It is very likely that something went wrong in the VASP relaxation calculation as the .out file is unreasonably short. If that is not the case, then modify the source code threshold constant.')
    
    return

def run_vasp_from_file(dirName, predefined_chgcar=None, run_type='ph'):
    dirName = checkPath(dirName)
    assert os.path.isdir(dirName)
    # TODO (function unnecessary for now)


