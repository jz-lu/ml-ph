from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_names import (
    CONFIG_BATCH_NAME, 
    START_BATCH_NAME, 
    CODE_DIR, MY_EMAIL, 
    SHIFT_NAME, 
    TOTAL_ENER_DIR_NAME, 
    POSCAR_NAME, 
    CONFIG_SUBDIR_NAME,
    ENERGIES
    )
from ___constants_compute import *
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from __class_input import InputData
from __class_Configuration import Configuration
from __directory_searchers import checkPath
from __dirModifications import mkdir
from __compute_properties import relax_solid
from __class_DataOutput import DataOutput
from multiprocessing import Pool
from random import randint
import os
import numpy as np
SINGLE_POOL = False

# Multiprocess for each displacement, since they are completely independent.
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
        with open(rtfile, 'w') as f:
            vdw = 'T' if user_input_settings.do_vdW else 'F'
            kpts = 'GAMMA' if user_input_settings.kpoints_is_gamma_centered else 'MP'
            randidx = randint(1,10000)
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -J %s\n'%(COMPUTE_JOBNAME))
            if USE_NODE_INDICATOR:
                f.write('#SBATCH -N %s\n'%(COMPUTE_NNODE))
            f.write('#SBATCH -n %s\n#SBATCH -t %s\n#SBATCH -p %s\n#SBATCH --mem-per-cpu=%s\n'%(COMPUTE_NCPU, COMPUTE_TIME, COMPUTE_PARTITIONS, COMPUTE_MEM_PER_CPU))
            f.write('#SBATCH -o shift_%A_%a.out\n#SBATCH -e shift_%A_%a.err\n')
            f.write('#SBATCH --mail-type=%s\n#SBATCH --mail-user=%s\n'%(COMPUTE_EMAIL_TYPE, COMPUTE_EMAIL_TO))
            f.write('source activate $HOME/%s\n'%(COMPUTE_ANACONDA_ENV))
            f.write('WDIR="%s${SLURM_ARRAY_TASK_ID}"\n'%(BASE_ROOT + CONFIG_SUBDIR_NAME))
            f.write('echo "WD: ${WDIR}"\n')
            f.write('ALLEGRO_DIR="%s"\n'%CODE_DIR)
            f.write('module list\nsource activate $HOME/%s\n'%(COMPUTE_ANACONDA_ENV))
            f.write('echo "RUNNING new array ridx = %d"\n'%randidx)
            f.write('python3 $ALLEGRO_DIR/start.py 0 $WDIR %s %s %s\n'%(vdw, kpts, ENERGIES))
            f.write('echo "FINISHED array ridx = %d"\n'%randidx)
        print('Main executable built.')

        for i, shpath in enumerate(base_root_subpaths):
            # with open(shpath + CONFIG_BATCH_NAME + str(i), 'w') as bshscr:
            #     bshscr.write('#!/bin/bash\n')
            #     bshscr.write('source activate $HOME/anaconda_env\n')
            #     bshscr.write('STR="%s"\n'%shpath)
            #     bshscr.write('echo "WD: ${STR}"\n')
            #     bshscr.write('ALLEGRO_DIR="%s"\n'%CODE_DIR)
            #     bshscr.write('python3 $ALLEGRO_DIR/start.py 0 $STR %s %s energies\n'%(vdw, kpts))
            bfile = shpath + SHIFT_NAME
            with open(bfile, 'w') as f:
                for val in configposcar_shift_tuple[i][0]:
                    f.write(str(val) + '\n')
            configposcar_shift_tuple[i][1].write_file(shpath + POSCAR_NAME) # write POSCAR to file
        print('All configurations executables built.')
        runcmd = 'sbatch --array=0-%d'%(nshifts-1) + ' ' + rtfile
        print('Running %s...'%runcmd)
        stream = os.popen(runcmd)
        print(stream.read())
        return None

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