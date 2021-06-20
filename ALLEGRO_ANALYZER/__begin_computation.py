from ____debug import DEBUGGING, DEBUG_NOTICE_MSG
from ___constants_config import GRID_SAMPLE_LOW, GRID_SAMPLE_HIGH
from ___constants_misc import NUM_AVAILABLE_CORES
from ___constants_vasp import Z_LAYER_SEP
from ___constants_compute import DEFAULT_NUM_LAYERS, MONOLAYER_JOBNAME, CONFIG_JOBNAME
from ___constants_names import (
    CONFIG_DATA_DIR, COB_NPY_NAME, 
    TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC, 
    POSCAR_NAME, POSCAR_CONFIG_NAMEPRE,
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME,
    PH
)
from __class_input import InputData
from __class_Configuration import Configuration
from __directory_searchers import checkPath, findFilesInDir
from __dirModifications import mkdir, copy, move
from __build_shscript import compute_configs
from __compute_properties import relax_solid, solve_electronic
from __class_DataOutput import DataOutput
from pymatgen.io.vasp.inputs import Poscar
from make_execscr import build_bash_exe
import os
import numpy as np

# Build configuration sampling from user input.
def begin_computation(user_input_settings):
    flg = user_input_settings.get_type_flag()
    BASE_ROOT = user_input_settings.get_base_root_dir()
    vdw = 'T' if user_input_settings.do_vdW else 'F'
    kpts = 'GAMMA' if user_input_settings.kpoints_is_gamma_centered() else 'MP'
    if flg == TYPE_RELAX_BASIC:
        print('Set to run standard single computation. Results to be stored to base root directory.')
        relax_solid(user_input_settings)
        return None # Indicate to start.py that there is no further data parsing necessary
    elif flg == TYPE_RELAX_CONFIG:
        print('Set to run parallel computations over grid sample in configuration space, defaulted to layer 1 (z = 0). Starting...')
        # Sample the grid here! Then pool them over to relax and run an iterator to get energy pairs for graphing
        # Num processes (in pool arg below) = number of grid points, i.e. (a, b, c) |-> a * b * c
        grid_size = GRID_SAMPLE_LOW # TODO: change LOW to HIGH in var name for super-accurate calculations
        init_interlayer_spacing = Z_LAYER_SEP

        sampling_set = Configuration.sample_grid(grid=grid_size)
        config = Configuration(BASE_ROOT)

        # Get a set of tuples (shift vector, poscar object) for each shift vector.
        configposcar_shift_tuple = config.build_config_poscar_set(sampling_set, init_interlayer_spacing)

        # Multiprocess over the configurations
        bze_points = compute_configs(BASE_ROOT, user_input_settings, configposcar_shift_tuple)

        # Parse the output of the configuration analysis.
        # Move home directory to the selected one.
        data_dir = user_input_settings.get_base_root_dir() + checkPath(CONFIG_DATA_DIR)
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        os.chdir(data_dir)
        print("Moved CWD to " + data_dir)

        # Get the 2D (without z vector) lattice basis from the POSCAR of the fixed layer,
        # then convert it into a change-of-basis matrix.
        lattice_basis = config.get_fixed_layer_poscar().as_dict()['structure']['lattice']['matrix'][:-1]
        cob_matrix = np.transpose([i[:-1] for i in lattice_basis])
        np.save(data_dir + COB_NPY_NAME, cob_matrix)
        print("Saved COB matrix for analysis at " + data_dir)

        # np.save(data_dir + 'bze', np.array(bze_points))
        # out = DataOutput(data_dir, bze_points, cob_matrix)
        # out.output_all_analysis()
        print("Configuration analysis (raw data file dump) complete. Returning data to `start.py`...")
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
        for i, p in enumerate(poscars):
            i = str(i)
            intra_path = checkPath(BASE_ROOT + MONOLAYER_DIR_NAME + i)
            mkdir(MONOLAYER_DIR_NAME + i, BASE_ROOT)
            copy(p, BASE_ROOT, newPath=intra_path)
            move(p, BASE_ROOT, newPath=BASE_ROOT+CONFIG_DIR_NAME)
            exepath = build_bash_exe(calc_type=TYPE_RELAX_BASIC, calc_list=[PH], outdir=intra_path,
                   compute_jobname=MONOLAYER_JOBNAME+i, vdw=vdw, kpts=kpts)
            if not DEBUGGING:
                print("Submitting monolayer job for layer " + i + "...")
                stream = os.popen('sbatch ' + exepath)
                print(stream.read())
        exepath = build_bash_exe(calc_type=TYPE_RELAX_CONFIG, calc_list=[PH], outdir=inter_path,
                   compute_jobname=CONFIG_JOBNAME, vdw=vdw, kpts=kpts)
        if not DEBUGGING:
            print("Submitting configuration job...")
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

