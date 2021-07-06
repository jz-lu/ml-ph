from ____debug import DEBUGGING, DEBUG_NOTICE_MSG
from ____exit_with_error import exit_with_error
from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from ___constants_compute import DEFAULT_NUM_LAYERS, MONOLAYER_JOBNAME, CONFIG_JOBNAME
from ___constants_names import (
    CONFIG_DATA_DIR, COB_NPY_NAME, LIDXS_NPY_NAME, 
    TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC, 
    POSCAR_NAME, POSCAR_CONFIG_NAMEPRE,
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME, 
    PH
)
from __class_input import InputData
from __class_Configuration import Configuration
from __class_RelaxerAPI import RelaxerAPI
from __directory_searchers import checkPath, findFilesInDir
from __dirModifications import mkdir, copy, move
from __build_shscript import compute_configs
from __compute_properties import relax_solid, solve_electronic
from pymatgen.io.vasp.inputs import Poscar
from make_execscr import build_bash_exe
import os
import numpy as np
RUN_RELAXER = False

# Build configuration sampling from user input.
def begin_computation(user_input_settings):
    flg = user_input_settings.get_type_flag()
    BASE_ROOT = user_input_settings.get_base_root_dir()
    vdw = user_input_settings.do_vdW
    kpts = user_input_settings.kpoints_is_gamma_centered
    if flg == TYPE_RELAX_BASIC:
        print('Set to run standard single computation. Results to be stored to base root directory.')
        relax_solid(user_input_settings)
        return None # Indicate to start.py that there is no further data parsing necessary
    elif flg == TYPE_RELAX_CONFIG:
        print('Set to run parallel computations over grid sample in configuration space, defaulted to layer 1 (z = 0). Starting...')
        grid_size = user_input_settings.get_cfg_grid_sz()
        print(f"Configuration grid: {grid_size}")
        init_interlayer_spacing = Z_LAYER_SEP
        config = Configuration(BASE_ROOT)
        data_dir = user_input_settings.get_base_root_dir() + checkPath(CONFIG_DATA_DIR)
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)

        sampling_set = None
        if user_input_settings.sampling_is_diagonal():
            print("Sampling line points uniformly along diagonal...")
            assert grid_size % 3 == 0, f"Number of points must be a multiple of 3 to properly assess high-symmetry points"
            sampling_set = Configuration.sample_line(npts=grid_size, basis=config.get_diagonal_basis())
            print(f"Sampled {grid_size} points along line {config.get_diagonal_basis()}")
        else:
            print("Sampling grid points along unit cell...")
            sampling_set = np.array(Configuration.sample_grid(grid=grid_size))
            if user_input_settings.get_tw_angle() is not None:
                print(f"Supplied twist angle: {user_input_settings.get_tw_angle()}")
            if RUN_RELAXER and user_input_settings.get_tw_angle() is not None:
                print("But first...computing relaxation in Julia...")
                # TODO cartesian or direct? Must convert if Cartesian using inverse COB
                u = RelaxerAPI(user_input_settings.get_tw_angle(), grid_size, data_dir).get_u()
                sampling_set = np.stack((sampling_set[:,:2] + u, sampling_set[:,2]), axis=1)
                print("Final b set:\n", sampling_set)


        # Get a set of tuples (shift vector, poscar object) for each shift vector.
        configposcar_shift_tuple = config.build_config_poscar_set(sampling_set, init_interlayer_spacing)

        # Multiprocess over the configurations
        bze_points = compute_configs(BASE_ROOT, user_input_settings, configposcar_shift_tuple)

        # Parse the output of the configuration analysis.
        # Move home directory to the selected one.
        os.chdir(data_dir)
        print("Moved CWD to " + data_dir)

        # Get the 2D (without z vector) lattice basis from the POSCAR of the fixed layer,
        # then convert it into a change-of-basis matrix.
        lattice_basis = config.get_fixed_layer_poscar().as_dict()['structure']['lattice']['matrix'][:-1]
        cob_matrix = np.transpose([i[:-1] for i in lattice_basis])
        np.save(data_dir + COB_NPY_NAME, cob_matrix)
        print("Saved COB matrix for analysis at " + data_dir)
        layer_idxs = config.get_layer_idxs()
        np.save(data_dir + LIDXS_NPY_NAME, layer_idxs)
        print(f"Saved layer indices {layer_idxs} for analysis at " + data_dir)

        print("Configuration analysis complete. Returning data to `start.py`...")
        print("**IMPORTANT**: when computation on all shifts complete, run `config_analyze.py` to get data analysis")
        return bze_points
    elif flg == TYPE_TWISTED_CONFIG:
        print('Splitting into 2 basic monolayer relaxations and 1 configuration calculation')
        poscars = sorted(findFilesInDir(BASE_ROOT, POSCAR_CONFIG_NAMEPRE, searchType='start'))
        assert len(poscars) > 1, "Must give at least 2 POSCARs, 1 per layer, for twist calculations"
        assert len(poscars) == DEFAULT_NUM_LAYERS, "Twist calculations for more than 2 layers not supported (yet)"

        print("Making new subdirectories, building I/O for monolayer and configuration calculations...")
        inter_path = checkPath(BASE_ROOT + CONFIG_DIR_NAME)
        mkdir(CONFIG_DIR_NAME, BASE_ROOT)
        clist = user_input_settings.get_calculation_list()
        for i, p in enumerate(poscars):
            i = str(i+1)
            intra_path = checkPath(BASE_ROOT + MONOLAYER_DIR_NAME + i)
            mkdir(MONOLAYER_DIR_NAME + i, BASE_ROOT)
            copy(p, BASE_ROOT, newPath=intra_path, newName=POSCAR_NAME)
            copy(p, BASE_ROOT, newPath=BASE_ROOT+CONFIG_DIR_NAME)
            exepath = build_bash_exe(calc_type='basic', calc_list=clist, outdir=intra_path,
                   compute_jobname=MONOLAYER_JOBNAME+i, vdw=vdw, kpts=kpts)
            if not DEBUGGING:
                os.chdir(intra_path)
                print("Submitting monolayer job for layer " + i + "...")
                stream = os.popen('sbatch ' + exepath)
                print(stream.read())
            else:
                print(DEBUG_NOTICE_MSG)
        exepath = build_bash_exe(calc_type='config', calc_list=clist, outdir=inter_path,
                   compute_jobname=CONFIG_JOBNAME, vdw=vdw, kpts=kpts, wdir=inter_path+CONFIG_SUBDIR_NAME,
                   compute_time='01:00:00', compute_ncpu='1', twist=user_input_settings.get_tw_angle(), 
                   sampling=user_input_settings.sampling) # this just kicks off a bunch of jobs, so it doesn't need any time
        if not DEBUGGING:
            os.chdir(inter_path)
            print(f"Submitting configuration job (wdir={inter_path+CONFIG_SUBDIR_NAME})...")
            stream = os.popen('sbatch ' + exepath)
            print(stream.read())
        else:
            print(DEBUG_NOTICE_MSG)
    else:
        print('Set to run analytical calculations directly, without relaxation')
        os.chdir(BASE_ROOT)
        poscars = findFilesInDir(BASE_ROOT, POSCAR_NAME, searchType='start')
        assert len(poscars) == 1
        p = Poscar.from_file(BASE_ROOT + poscars[0])
        solve_electronic(user_input_settings, poscar=p)
        return None

