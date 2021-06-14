from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_names import CONFIG_BATCH_NAME, START_BATCH_NAME, CODE_DIR, MY_EMAIL, SHIFT_NAME, TOTAL_ENER_DIR_NAME
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from __class_input import InputData
from __class_Configuration import Configuration
from __directory_searchers import checkPath
from __dirModifications import mkdir
from __compute_properties import relax_solid
from __class_DataOutput import DataOutput
from multiprocessing import Pool
import os
import numpy as np
SINGLE_POOL = False

# Multiprocess for each displacement, since they are completely independent.
def compute_configs(BASE_ROOT, user_input_settings, configposcar_shift_tuple):
    BASE_ROOT = checkPath(BASE_ROOT)
    SUBDIRNAMES = 'shift_'
    base_root_subpaths = []
    print('Creating subdirectories to store VASP calculations for each configuration...')
    for i in range(len(configposcar_shift_tuple)):
        new_subdir_name = '%s%d/'%(SUBDIRNAMES, i)
        mkdir(new_subdir_name, BASE_ROOT)
        base_root_subpaths.append(BASE_ROOT + new_subdir_name)
    
    if SINGLE_POOL: # Run the whole job (all configs) in one go now
        return one_job_run(user_input_settings, configposcar_shift_tuple, base_root_subpaths)
    else: # Make a new job for each configuration and submit bash scripts to the cluster
        rtfile = BASE_ROOT + START_BATCH_NAME
        print('Building configuration space job array executables...')
        with open(rtfile, 'w') as f:
            f.write('#!/bin/bash\n')
            f.write('#SBATCH -J shifts\n')
            f.write('#SBATCH -n 1\n#SBATCH -t 8:00:00\n#SBATCH -p kaxiras,shared\n#SBATCH --mem-per-cpu=5000\n')
            f.write('#SBATCH -o shift_%A_%a.out\n#SBATCH -e shift_%A_%a.err\n')
            f.write('#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=%s\n'%MY_EMAIL)
            f.write('module list\nsource activate $HOME/anaconda_env\n')
            
            f.write('echo "RUNNING array idx %A_%a"\n')
            f.write(START_BATCH_NAME + ' ' + BASE_ROOT + SUBDIRNAMES + '"${SLURM_ARRAY_TASK_ID}"' + CONFIG_BATCH_NAME + '"${SLURM_ARRAY_TASK_ID}"')
            f.write('echo "FINISHED array idx %A_%a"\n')
        print('Main executable built.')

        vdw = 'T' if user_input_settings.do_vdW else 'F'
        kpts = 'GAMMA' if user_input_settings.kpoints_is_gamma_centered else 'MP'
        for i, shpath in enumerate(base_root_subpaths):
            with open(shpath + CONFIG_BATCH_NAME + str(i), 'w') as bshscr:
                bshscr.write('#!/bin/bash\n')
                bshscr.write('source activate $HOME/anaconda_env\n')
                bshscr.write('STR="%s"\n'%BASE_ROOT)
                bshscr.write('echo "WD: ${STR}"\n')
                bshscr.write('ALLEGRO_DIR="%s"\n'%CODE_DIR)
                bshscr.write('python3 $ALLEGRO_DIR/start.py 0 $STR %s %s energies\n'%(vdw, kpts))
            bfile = shpath + SHIFT_NAME
            with open(bfile, 'w') as f:
                for val in configposcar_shift_tuple[i][0]:
                    f.write(str(val))
        print('All configurations executables built.')
        print('Running executables...')
        stream = os.popen('sbatch ' + rtfile)
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