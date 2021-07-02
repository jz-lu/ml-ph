from ____debug import DEBUGGING, DEBUG_NOTICE_MSG
from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_names import (
    CONFIG_BATCH_NAME, 
    START_BATCH_NAME, 
    CODE_DIR, MY_EMAIL, 
    SHIFT_NAME, 
    TOTAL_ENER_DIR_NAME, 
    POSCAR_NAME, 
    CONFIG_SUBDIR_NAME,
    ENERGIES,
    PHDISP_STATIC_NAME,
    TYPE_RELAX_BASIC
    )
from ___constants_compute import *
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from __class_input import InputData
from __class_Configuration import Configuration
from __directory_searchers import checkPath
from __dirModifications import mkdir
from __compute_properties import relax_solid
from multiprocessing import Pool
from make_execscr import build_bash_exe
import os
import numpy as np
SINGLE_POOL = False

# Compute configurations by submitting an array of independent jobs to the computing cluster
def compute_configs(BASE_ROOT, user_input_settings, configposcar_shift_tuple):
    BASE_ROOT = checkPath(BASE_ROOT)
    base_root_subpaths = []
    nshifts = len(configposcar_shift_tuple)
    print('Creating subdirectories to store VASP calculations for each configuration...')
    for i in range(nshifts):
        new_subdir_name = '%s%d/'%(CONFIG_SUBDIR_NAME, i)
        mkdir(new_subdir_name, BASE_ROOT)
        base_root_subpaths.append(BASE_ROOT + new_subdir_name)
    
    if SINGLE_POOL: # Run the whole job (all configs) in one go now
        return one_job_run(user_input_settings, configposcar_shift_tuple, base_root_subpaths)
    else: # Make a new job for each configuration and submit bash scripts to the cluster
        rtfile = BASE_ROOT + START_BATCH_NAME
        print('Building configuration space job array executables...')
        vdw = 'T' if user_input_settings.do_vdW else 'F'
        kpts = 'GAMMA' if user_input_settings.kpoints_is_gamma_centered else 'MP'
        compute_time = '16:00:00' if vdw == 'T' else '08:00:00'
        compute_ncpu = '16' if vdw == 'T' else '8'
        build_bash_exe(calc_type=TYPE_RELAX_BASIC, outdir=BASE_ROOT, 
                        compute_jobname=DIAG_JOBNAME if user_input_settings.sampling_is_diagonal() else COMPUTE_JOBNAME, 
                        compute_time=compute_time, vdw=vdw, kpts=kpts, fname=START_BATCH_NAME, as_arr=True, compute_ncpu=compute_ncpu, 
                        wdir=BASE_ROOT+CONFIG_SUBDIR_NAME, calc_list=user_input_settings.get_raw_calculation_list())

        for i, shpath in enumerate(base_root_subpaths):
            bfile = shpath + SHIFT_NAME
            with open(bfile, 'w') as f:
                for val in configposcar_shift_tuple[i][0]:
                    f.write(str(val) + '\n')
            configposcar_shift_tuple[i][1].write_file(shpath + POSCAR_NAME) # write POSCAR to file
        print('Configurations written to file.')
        if not DEBUGGING:
            runcmd = 'sbatch --array=0-%d'%(nshifts-1) + ' ' + rtfile
            print('Running %s...'%runcmd)
            stream = os.popen(runcmd)
            print(stream.read())
        else:
            print(DEBUG_NOTICE_MSG)
        return None

# Multiprocess directly for each displacement, since they are completely independent.
def one_job_run(user_input_settings, configposcar_shift_tuple, base_root_subpaths):
    # Create a persistent pool of processors computing the configuration forces.
    print('Constructing multiprocessing worker pool...')
    pool = Pool(NUM_AVAILABLE_CORES) # Parallelize as much as possible.
    print('Pool loaded.')

    # Run asynchronous isolated calculations for each displacement over the pool.
    print('Starting parallel computation of configurations...')
    bze_points = [pool.apply_async(relax_solid, args=(user_input_settings, 
                                                      configposcar_shift_tuple[i][1], 
                                                      configposcar_shift_tuple[i][0], 
                                                      base_root_subpaths[i])).get() for i in range(len(configposcar_shift_tuple))]

    return bze_points


