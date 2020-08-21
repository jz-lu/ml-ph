import sys, os, copy
from pymatgen.io.vasp.inputs import Kpoints, Poscar, VaspInput
from pymatgen.io.vasp.outputs import Outcar

from ____exit_with_error import exit_with_error

from ___constants_misc import *
from ___constants_names import *
from ___constants_vasp import VASP_OUTFILE_LEN_THRESHOLD

from __directory_searchers import checkPath
from __dirModifications import mkdir
from __build_inputs import buildInitialInputs
from __run_vasp import run_vasp
from __get_energy_analysis import get_energies
from __postprocess_relaxation import postProcess_relaxation
from __class_input import InputData
from __class_CarCollector import CarCollector

async def relax_solid(user_input_settings):
    # Collect input files
    vasp_input_initializers = CarCollector(user_input_settings.ROOT, user_input_settings.do_vdW, user_input_settings.kpoints_is_gamma_centered, user_input_settings.need_line_kpoints())
    init_vasp_obj = vasp_input_initializers.build_relaxation_input(print_all_info=True)

    # Before we build the relaxation calculation inputs below, we need a sanity check as to whether the calculations can be done
    # The only thing we need to check besides poscar validity (code will handle other input files if none are given
    # is the existence of a line kpoints file if we want an electronic band structure calculation
    kpoints_line = CarCollector.get_line_kpoints(user_input_settings.ROOT, user_input_settings.need_line_kpoints())

    # We want to run the relaxation on a separate subfolder, so everything is organized
    mkdir(RELAXATION_DIR_NAME, user_input_settings.ROOT)
    DIR_RELAXATION = checkPath(user_input_settings.ROOT + RELAXATION_DIR_NAME)
    print('Created new folder %s to store relaxation calculations.'%(DIR_RELAXATION))

    # Similarly we want to run the analyses post-relaxation in a separate subfolder
    mkdir(ANALYSIS_DIR_NAME, user_input_settings.ROOT)
    DIR_ANALYSIS = checkPath(user_input_settings.ROOT + ANALYSIS_DIR_NAME)
    print('Created new folder %s to store analysis calculations.'%(DIR_ANALYSIS))

    # Call the relaxation
    print('Running VASP relaxation calculations...results to be sent to %s'%(DIR_RELAXATION))
    # print(init_vasp_obj)
    run_vasp(user_input_settings.ROOT, DIR_RELAXATION)
    print('VASP relaxation calculations complete.')

    # Begin post-processing the relaxation calculations, handled by the postprocessing module
    # Step 1: do the total energy
    if user_input_settings.do_energy_calculation():
        energy_pair = get_energies(DIR_RELAXATION, DIR_ANALYSIS, writeOut=True)
        
    if not user_input_settings.do_nonenergy_calculations():
        print('Successfully completed total energy calculation, which was the only specified calculation. Exiting...')
        sys.exit()

    # Step 2, for ele/ph calculations, we transfer to postprocessing module
    postProcess_relaxation(DIR_ANALYSIS, DIR_RELAXATION, init_vasp_obj, user_input_settings.get_calculation_list(), kpoints_line)
    return energy_pair # if we ever need anything to return up to config.py from this it is the energies