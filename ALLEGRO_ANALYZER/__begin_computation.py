from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from __class_input import InputData
from __class_Configuration import Configuration
from __directory_searchers import checkPath
from __dirModifications import mkdir
from multiprocessing import Pool
from __compute_properties import relax_solid

# Multiprocess for each displacement, since they are completely independent.
def branch_and_compute(BASE_ROOT, user_input_settings, configposcar_shift_tuple):
    BASE_ROOT = checkPath(BASE_ROOT)
    base_root_subpaths = []
    print('Creating subdirectories to store VASP calculations for each configuration...')
    for i in range(len(configposcar_shift_tuple)):
        new_subdir_name = 'shift_%d'%(i+1)
        mkdir(new_subdir_name, BASE_ROOT)
        base_root_subpaths.append(BASE_ROOT + new_subdir_name)

    # Create a persistent pool of processors computing the configuration forces.
    pool = Pool(NUM_AVAILABLE_CORES) # Parallelize as much as possible.

    # Run asynchronous isolated calculations for each displacement over the pool.
    print('Starting parallel computation of configurations.')
    bze_points = [pool.apply_async(relax_solid, args=(user_input_settings, 
                                                      configposcar_shift_tuple[i][1], 
                                                      configposcar_shift_tuple[i][0], 
                                                      base_root_subpaths[i])).get() for i in range(len(configposcar_shift_tuple))]

    return bze_points

# Build configuration sampling from user input.
def begin_computation(user_input_settings):
    if user_input_settings.get_type_flag() == 0:
        print('Set to run standard single computation. Results to be stored to base root directory.')
        relax_solid(user_input_settings)
        return None
    else:
        print('Set to run parallel computations over grid sample in configuration space, defaulted to layer 1 (z = 0). Starting...')
        # Sample the grid here! Then pool them over to relax and run an iterator to get energy pairs for graphing
        # Num processes (in pool arg below) = number of grid points, i.e. (a, b, c) |-> a * b * c
        grid_size = GRID_SAMPLE_LOW # NOTE: change to ..._HIGH for super-accurate calculations
        BASE_ROOT = user_input_settings.get_base_root_dir()
        init_interlayer_spacing = Z_LAYER_SEP

        sampling_set = Configuration.sample_grid(grid=grid_size)
        config = Configuration(BASE_ROOT)

        configposcar_shift_tuple = config.build_config_poscar_set(sampling_set, init_interlayer_spacing)

        # Multiprocess over the configurations
        bze_points = branch_and_compute(BASE_ROOT, user_input_settings, configposcar_shift_tuple)
        return bze_points
    
