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
        vdw = 'T' if user_input_settings.do_vdW in ['T', True] else 'F'
        kpts = 'GAMMA' if user_input_settings.kpoints_is_gamma_centered else 'MP'
        compute_time = '12:00:00' if vdw == 'T' else '16:00:00'
        compute_ncpu = '48' if vdw == 'T' else '24' # change as necessary
        compute_jobname = DIAG_JOBNAME if user_input_settings.sampling_is_diagonal() else COMPUTE_JOBNAME
        compute_jobname = user_input_settings.passname() + compute_jobname
        compute_jobname = '-' + compute_jobname
        if user_input_settings.run_relaxer:
            compute_jobname = 'r' + compute_jobname
        if vdw == 'T':
            compute_jobname = 'v' + compute_jobname
        if user_input_settings.sampling_is_interlayer():
            compute_jobname = 'z' + compute_jobname
        elif user_input_settings.sampling == 'low':
            compute_jobname = 'l' + compute_jobname
        else:
            compute_jobname = 'h' + compute_jobname
        no_ionic_step = False
        if user_input_settings.sampling_is_interlayer():
            no_ionic_step = True
        build_bash_exe(calc_type='basic', outdir=BASE_ROOT,
                        compute_jobname=compute_jobname, compute_time=compute_time, vdw=vdw, kpts=kpts, fname=START_BATCH_NAME, 
                        as_arr=True, compute_ncpu=compute_ncpu, wdir=BASE_ROOT+CONFIG_SUBDIR_NAME, 
                        calc_list=user_input_settings.get_raw_calculation_list(), 
                        compute_partitions=COMPUTE_PARTITIONS, 
                        passname=user_input_settings.passname(), pass_idx=True, 
                        fcut=user_input_settings.fcut, ediff0=user_input_settings.ediff0, 
                        super_dim=user_input_settings.get_super_dim(), 
                        no_ionic_step=no_ionic_step)

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


