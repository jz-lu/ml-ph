from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from __class_input import InputData
from __class_Configuration import Configuration
from __directory_searchers import checkPath
from __dirModifications import mkdir
from multiprocessing import Pool
from __compute_properties import relax_solid
from __class_DataOutput import DataOutput
import os
import numpy as np
from sys import exit

# Multiprocess for each displacement, since they are completely independent.
def branch_and_compute(BASE_ROOT, user_input_settings, configposcar_shift_tuple):
    BASE_ROOT = checkPath(BASE_ROOT)
    base_root_subpaths = []
    print('Creating subdirectories to store VASP calculations for each configuration...')
    for i in configposcar_shift_tuple:
        new_subdir_name = 'shift_' + str(tuple(i[0]))
        mkdir(new_subdir_name, BASE_ROOT)
        base_root_subpaths.append(BASE_ROOT + new_subdir_name)

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

# Build configuration sampling from user input.
def begin_computation(user_input_settings):
    if user_input_settings.get_type_flag() == 0:
        print('Set to run standard single computation. Results to be stored to base root directory.')
        relax_solid(user_input_settings)
        return None # Indicate to start.py that there is no further data parsing necessary.
    else:
        print('Set to run parallel computations over grid sample in configuration space, defaulted to layer 1 (z = 0). Starting...')
        # Sample the grid here! Then pool them over to relax and run an iterator to get energy pairs for graphing
        # Num processes (in pool arg below) = number of grid points, i.e. (a, b, c) |-> a * b * c
        grid_size = GRID_SAMPLE_LOW # NOTE: change LOW to HIGH oin var name for super-accurate calculations
        BASE_ROOT = user_input_settings.get_base_root_dir()
        init_interlayer_spacing = Z_LAYER_SEP

        sampling_set = Configuration.sample_grid(grid=grid_size)
        config = Configuration(BASE_ROOT)

        # Get a set of tuples (shift vector, poscar object) for each shift vector.
        configposcar_shift_tuple = config.build_config_poscar_set(sampling_set, init_interlayer_spacing)

        # Multiprocess over the configurations
        bze_points = branch_and_compute(BASE_ROOT, user_input_settings, configposcar_shift_tuple)

        # Parse the output of the configuration analysis.
        # Move home directory to the selected one.
        data_dir = user_input_settings.get_base_root_dir() + "raw_data/"
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        os.chdir(data_dir)
        print("Moved CWD to " + data_dir)

        # Get the 2D (without z vector) lattice basis from the POSCAR of the fixed layer,
        # then convert it into a change-of-basis matrix.
        lattice_basis = config.get_fixed_layer_poscar().as_dict()['structure']['lattice']['matrix'][:-1]
        cob_matrix = np.transpose([i[:-1] for i in lattice_basis])

        out = DataOutput(data_dir, bze_points, cob_matrix)
        out.output_all_analysis()
        print("Configuration analysis (raw table export, plotting) complete. Returning data to start.py...")

        return bze_points
    
