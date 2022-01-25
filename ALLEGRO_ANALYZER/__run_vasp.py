# Main file to be executed by user
import os
import numpy as np
import shutil
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
    KPOINTS_NAME, 
    CODE_DIR
)
from ___constants_misc import ERR_VASP_RUN_RELAX, ERR_VASP_NOT_CONVERGED
from ___constants_vasp import (
    VASP_OUTFILE_LEN_THRESHOLD, VASP_MAX_CONVERGENCE_ATTEMPTS, 
    IBRION
)
from __class_CarCollector import CarCollector
from __directory_searchers import checkPath
from __check_convergence import check_if_not_converged

# Implement the loosening part of the loosen-tighten algorithm to 
# maximize efficient convergence at desired EDIFF in a relaxation.
def run_vasp_relaxation(vaspObj, dirName, outfile_name, errfile_name):
    loopCounter = 1 # Counts number of times we try to relax. We stop at VASP_MAX_CONVERGENCE_ATTEMPTS

    # Run vasp and check convergence
    vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=outfile_name, err_file=errfile_name)
    not_converged = check_if_not_converged(dirName, outfile_name, max_nsw=vaspObj['INCAR']['NSW'])
    if not not_converged:
        print(f"LT: step at {vaspObj['INCAR']['EDIFF']} complete")

    while not_converged:
        print('Relaxation has not converged. Proceeding to rerun relaxation with updated ionic positions...')
        print('Loop counter: ' + str(loopCounter))
        loopCounter += 1
        if loopCounter > VASP_MAX_CONVERGENCE_ATTEMPTS['L']:
            if vaspObj['INCAR']['EDIFFG'] < 1e-3:
                vaspObj['INCAR']['EDIFFG'] = round(10 * vaspObj['INCAR']['EDIFFG'], 12)
                vaspObj['INCAR']['EDIFF'] = round(10 * vaspObj['INCAR']['EDIFF'], 12)
                print(f"LT: loosening EDIFF to {vaspObj['INCAR']['EDIFF']}")
                print("Loop counter reset.")
                loopCounter = 1
            else:
                print("Unable to further loosen EDIFFG threshold.")
                print(ERR_VASP_NOT_CONVERGED)
                break
        
        # Update POSCAR to be the new CONTCAR
        new_poscar = Poscar.from_file(dirName + CONTCAR_NAME)
        vaspObj['POSCAR'] = new_poscar 
        print('CONTCAR copied to POSCAR. Rerunning relaxation...')
        # No need to delete anything, the files will be overwritten. Just start the calculation.
        vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=outfile_name, err_file=errfile_name)
        print('Run complete.')
        not_converged = check_if_not_converged(dirName, outfile_name, max_nsw=vaspObj['INCAR']['NSW'])
        if not not_converged:
            print(f"LT: step at {vaspObj['INCAR']['EDIFF']} complete")
    return vaspObj

def run_vasp(vaspObj, dirName, predefined_chgcar=None, run_type='relax', edinit=None, force_cut=False):
    assert (edinit is not None) or run_type != 'relax', f"Initial EDIFF not given for relaxation"
    if DEBUGGING:
        print(DEBUG_NOTICE_MSG)
        return
    dirName = checkPath(dirName)
    
    # Copy vdW kernel into the right folder
    try:
        shutil.copyfile(CODE_DIR + 'vdw_kernel.bindat', dirName)
    except BaseException as e:
        print(e)
    
    # pymatgen doesn't handle other files well. We will manually print chgcar in the outdir so VASP will take it when it runs.
    if predefined_chgcar != None:
        predefined_chgcar.write_file(dirName + CHGCAR_NAME)
    
    outfile_name = VASP_RUN_OUT_NAME%(run_type)
    errfile_name = VASP_RUN_ERR_NAME%(run_type)
    if run_type == 'dfpt':
        print("Using DFPT algorithm", flush=True)
        vaspObj['INCAR']['IBRION'] = IBRION['dfpt']
        vaspObj['INCAR']['NSW'] = 1
        if 'NPAR' in list(vaspObj['INCAR'].keys()):
            vaspObj['INCAR'].pop('NPAR')

    # No relaxation -> single run, no need to check for ionic convergence
    if run_type != 'relax':
        vaspObj.run_vasp(run_dir=dirName, vasp_cmd=[BATCH_FILE_PATH], output_file=outfile_name, err_file=errfile_name)

    # Relaxation -> need to make sure it converges
    else:
        EDIFF0 = vaspObj['INCAR']['EDIFF']
        vaspObj['INCAR']['NPAR'] = 4
        sgn = -1 if force_cut else 1
        vaspObj['INCAR']['EDIFF'] = edinit
        vaspObj['INCAR']['EDIFFG'] = vaspObj['INCAR']['EDIFFG'] = sgn * edinit # round(10 * edinit, 12)
        print(f"Starting initial LT-algorithm run at EDIFF={edinit}", flush=True)
        vaspObj = run_vasp_relaxation(vaspObj, dirName, outfile_name, errfile_name)
        LT_converged = False
        tight_loopctr = 1
        print("LT: switching to quasi-Newton method descent")
        vaspObj['INCAR']['IBRION'] = 1
        # vaspObj['INCAR']['POTIM'] = 0.1
        while tight_loopctr <= VASP_MAX_CONVERGENCE_ATTEMPTS['T']:
            EDIFF_here = vaspObj['INCAR']['EDIFF']
            if np.isclose(EDIFF_here, EDIFF0) or EDIFF_here <= EDIFF0 or round(EDIFF_here, 12) <= 1e-8:
                LT_converged = True
                print(f"[{tight_loopctr}] LT: relaxation has fully converged")
                break
            else:
                print(f"[{tight_loopctr}] LT: EDIFF={EDIFF_here} not to thres level {EDIFF0} yet")
            vaspObj['INCAR']['EDIFFG'] = round(vaspObj['INCAR']['EDIFFG'] / 10, 12)
            vaspObj['INCAR']['EDIFF'] = round(vaspObj['INCAR']['EDIFF'] / 10, 12)
            print(f"LT: tightening EDIFF by one order to {vaspObj['INCAR']['EDIFF']}")
            new_poscar = Poscar.from_file(dirName + CONTCAR_NAME)
            vaspObj['POSCAR'] = new_poscar
            vaspObj = run_vasp_relaxation(vaspObj, dirName, outfile_name, errfile_name)

            if round(vaspObj['INCAR']['EDIFF'], 12) > EDIFF_here:
                tight_loopctr += 1
        if not LT_converged:
            print(f"LT: too many runs failed with {tight_loopctr} T-runs. Keeping best EDIFF={vaspObj['INCAR']['EDIFF']}")
        
    # Check the most obvious error in the run automatically
    if os.stat(dirName + outfile_name).st_size < VASP_OUTFILE_LEN_THRESHOLD:
        exit_with_error('It is very likely that something went wrong in the VASP relaxation calculation as the .out file is unreasonably short. If that is not the case, then modify the source code threshold constant.')
    
    return

def run_vasp_from_file(dirName, predefined_chgcar=None, run_type='ph'):
    dirName = checkPath(dirName)
    assert os.path.isdir(dirName)
    # TODO (function unnecessary for now)


