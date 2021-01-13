from ____exit_with_error import exit_with_error
from ___constants_names import RELAXATION_DIR_NAME, ANALYSIS_DIR_NAME, CONTCAR_NAME
from __directory_searchers import checkPath
from __dirModifications import mkdir
from __build_inputs import buildInitialInputs
from __run_vasp import run_vasp
from __get_energy_analysis import get_energies
from __postprocess_relaxation import postProcess_relaxation
from __class_CarCollector import CarCollector
import sys, os, copy
from pymatgen.io.vasp.inputs import Kpoints, Poscar, VaspInput
from pymatgen.io.vasp.outputs import Outcar

def relax_solid(user_input_settings, poscar=None, shift=None, user_inputted_root=None):
    # Get the right directory to import the relevant files.
    ROOT = user_input_settings.get_base_root_dir()
    if user_inputted_root != None:
        ROOT = user_inputted_root
    ROOT = checkPath(ROOT)

    # Collect input files
    vasp_input_initializers = CarCollector(ROOT, user_input_settings.do_vdW, user_input_settings.kpoints_is_gamma_centered, user_input_settings.need_line_kpoints(), poscar=poscar)
    init_vasp_obj = vasp_input_initializers.build_relaxation_input(print_all_info=True)

    # Before we build the relaxation calculation inputs below, we need a sanity check as to whether the calculations can be done
    # The only thing we need to check besides poscar validity (code will handle other input files if none are given
    # is the existence of a line kpoints file if we want an electronic band structure calculation
    kpoints_line = CarCollector.get_line_kpoints(ROOT, user_input_settings.need_line_kpoints())

    output_ROOT = 

    # We want to run the relaxation on a separate subfolder, so everything is organized
    mkdir(RELAXATION_DIR_NAME, ROOT)
    DIR_RELAXATION = checkPath(ROOT + RELAXATION_DIR_NAME)
    print('Created new folder %s to store relaxation calculations.'%(DIR_RELAXATION))

    # Similarly we want to run the analyses post-relaxation in a separate subfolder
    mkdir(ANALYSIS_DIR_NAME, ROOT)
    DIR_ANALYSIS = checkPath(ROOT + ANALYSIS_DIR_NAME)
    print('Created new folder %s to store analysis calculations.'%(DIR_ANALYSIS))

    # Call the relaxation
    print('Running VASP relaxation calculations...results to be sent to %s'%(DIR_RELAXATION))
    # print(init_vasp_obj)
    run_vasp(ROOT, DIR_RELAXATION)
    print('VASP relaxation calculations complete.')

    # Begin post-processing the relaxation calculations, handled by the postprocessing module
    # Step 1: do the total energy
    if user_input_settings.do_energy_calculation():
        energy_pair = get_energies(DIR_RELAXATION, DIR_ANALYSIS, writeOut=True)
        
    if user_input_settings.do_nonenergy_calculations():
        # Step 2, for ele/ph calculations, we transfer to postprocessing module
        postProcess_relaxation(DIR_ANALYSIS, DIR_RELAXATION, init_vasp_obj, user_input_settings.get_calculation_list(), kpoints_line)
    else:
        print('Successfully completed total energy calculation, which was the only specified calculation. Exiting relaxation...')

    if shift == None:
        sys.exit()
    else:
        bze_point = (shift, CarCollector.get_interlayer_spacing(DIR_RELAXATION), energy_pair[0])
        return bze_point

