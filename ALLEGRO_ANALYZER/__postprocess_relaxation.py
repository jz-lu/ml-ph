import os
import subprocess
import matplotlib.pyplot as plt # pylint: disable=import-error
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput # pylint: disable=import-error 
from pymatgen.io.vasp.outputs import Chgcar # pylint: disable=import-error
from pymatgen.core.structure import Structure # pylint: disable=import-error
from __build_inputs import buildPotcar
from __input_modifiers import modifyIncar, modifyKpointsMesh, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy modules to the command line pipeline for basic linux commands
import copy
from __directory_searchers import checkPath
from __run_vasp import run_vasp
from __ph_processing import ph_prepare_for_analysis, ph_preprocess
from __get_ph_analysis import ph_get_dos, ph_get_band
from __get_eledos_analysis import get_eledos_analysis
from __get_eleband_analysis import get_eleband_analysis
from __get_elecombined_analysis import get_elecombined_analysis
from __cleanup import cleanRelevantFiles
from ___constants_vasp import NEDOS, ICHARG, SIGMA, PHONOPY_GRID_DENSITY, PHONOPY_GRID_SHIFT
from ___constants_names import (
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    POSCAR_NAME, CHGCAR_NAME, CONTCAR_NAME, 
    ELEDOS, ELEBAND, PH, PHDOS, PHBAND, ENERGIES, 
    OUTPUT_DIR_NAME, DOSCAR_NAME, THIS_DIR
)
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG
from ___constants_phonopy import POSCAR_UNIT_NAME
import os


# Process the command-line arguments
def postProcess_relaxation(outDirName, relaxation_dirName, unrelaxed_vaspObj, user_input_settings, kpoints_line): # vasp input object generated at beginning of relaxation
    print('Postprocessing VASP relaxation beginning...')
    outDirName = checkPath(outDirName)
    relaxation_dirName = checkPath(relaxation_dirName)

    # Get the relevant objects for the next set of calculations
    calculation_list = sorted(list(user_input_settings.get_calculation_list()))
    chgcar = None
    if ELEDOS in calculation_list or ELEBAND in calculation_list:
        chgcar = Chgcar.from_file(relaxation_dirName + CHGCAR_NAME) # Need to retrieve the last charge density to get good calculations
    poscar_relaxed = Poscar.from_file(relaxation_dirName + CONTCAR_NAME) # Get the relaxed positions from CONTCAR
    # Make deep copies here since we'll need the original
    relaxation_incar_1 = copy.deepcopy(unrelaxed_vaspObj['INCAR']) # For the inspective code reviewer: these are pymatgen's hard strings so no need to collect them
    relaxation_incar_2 = copy.deepcopy(relaxation_incar_1)
    potcar = buildPotcar(outDirName, poscar_relaxed) # outdirname is irrelevant here
    
    # Incar changes for various calculations
    incar_selfcon = getSelfConNoRelIncar(relaxation_incar_1)
    incar_nonselfcon = getNonSelfConNoRelIncar(relaxation_incar_2)
    print('[DEBUGMSG] incar selfcon:', incar_selfcon)
    print('[DEBUGMSG] incar non selfcon:', incar_nonselfcon)
    
    # Kpoints needs to be modified to have a denser sampling for nonrelaxation calculations -> more precision
    kpoints_mesh_nonrelax = copy.deepcopy(unrelaxed_vaspObj['KPOINTS'])
    kpoints_mesh_nonrelax = modifyKpointsMesh(kpoints_mesh_nonrelax)
    # For phonopy it might need to be a little denser but not necessarily as dense as dos
    kpoints_mesh_ph = copy.deepcopy(unrelaxed_vaspObj['KPOINTS'])
    kpoints_mesh_ph = modifyKpointsMesh(kpoints_mesh_ph, meshDensity=PHONOPY_GRID_DENSITY, totShift=PHONOPY_GRID_SHIFT)

    # Objects useful for plotting data
    eledos_obj = None # pylint: disable=unused-variable
    eleband_obj = None # pylint: disable=unused-variable
    combinedPlot = None # pylint: disable=unused-variable
    
    # Since all processing for phonopy is done at once, we should flag so that if we want phband and phdos we don't preprocess twice
    ph_has_preprocessed = False
    DIR_PHONOPY = None # This will be set during first processing

    # We want eleband to follow eledos if both are calculated since we can get a better CHGCAR for band from DOS
    # So we put it in the front of the array
    eledos_has_run = False
    if ELEDOS in calculation_list:
        calculation_list.insert(0, ELEDOS)
        calculation_list = list(dict.fromkeys(calculation_list)) # remove the duplicate, the first will be the one

    print('Ordered list of calculations:', calculation_list)

    # Parse command line and direct necessary function calls
    for i in calculation_list:
        if i == ELEDOS:
            print('Now running electronic DOS calculations.')
            mkdir(i, outDirName) # Create a subfolder for the analysis
            DIR_ELEDOS = checkPath(outDirName + ELEDOS)

            # os.chdir(DIR_ELEDOS) # Go to that working directory

            mkdir(OUTPUT_DIR_NAME, DIR_ELEDOS) # subsubfolder for result storage
            DIR_ELEDOS_RESULTS = checkPath(DIR_ELEDOS + OUTPUT_DIR_NAME)

            # In addition to self-consistent calculation, we need to also add on the NEDOS parameter to sample the space for eledos
            incar_eledos = copy.deepcopy(incar_selfcon)
            incar_eledos = modifyIncar(incar_eledos, addArr=[('NEDOS', NEDOS)])

            # We will need the standard VASP IO Object plus CHGCAR.
            eledos_vasp_obj = VaspInput(incar_eledos, kpoints_mesh_nonrelax, poscar_relaxed, potcar)
            print('[DEBUGMSG] eledos vasp incar: ', eledos_vasp_obj['INCAR'])

            print('Running VASP to get data for electronic DOS...')
            # Run vasp nonrelaxation, self-consistent
            run_vasp(eledos_vasp_obj, DIR_ELEDOS, predefined_chgcar=chgcar, run_type=ELEDOS)
            eledos_has_run = True

            # Get the analysis plots and data
            print('Parsing VASP run and retrieving DOS data...')
            eledos_obj = get_eledos_analysis(DIR_ELEDOS, DIR_ELEDOS_RESULTS, poscar_relaxed, extract_raw_data=True, extract_plot=True)
            if os.path.isfile(DIR_ELEDOS + DOSCAR_NAME):
                move(DOSCAR_NAME, DIR_ELEDOS, DIR_ELEDOS_RESULTS)
        elif i == ELEBAND:
            print('Now running electronic band structure calculations.')
            mkdir(i, outDirName) # Create a subfolder for the analysis
            DIR_ELEBAND = checkPath(outDirName + ELEBAND)

            # os.chdir(DIR_ELEBAND) # Go to that working directory

            mkdir(OUTPUT_DIR_NAME, DIR_ELEBAND) # subsubfolder for result storage
            DIR_ELEBAND_RESULTS = checkPath(DIR_ELEBAND + OUTPUT_DIR_NAME)

            # Update to best possible CHGCAR
            if eledos_has_run:
                print('Electronic DOS has already run. Fetching updated charge densities from dos run...')
                chgcar = Chgcar.from_file(checkPath(outDirName + ELEDOS) + CHGCAR_NAME)

            # We use the line kpoints file that we imported in the command line parsing start.py
            eleband_vasp_obj = VaspInput(incar_nonselfcon, kpoints_line, poscar_relaxed, potcar)
            print('Running VASP for electronic band calculations...')
            run_vasp(eleband_vasp_obj, DIR_ELEBAND, predefined_chgcar=chgcar, run_type=ELEBAND)

            # Get the analysis plots and data
            print('Parsing VASP run and retrieving band structure data...')
            eleband_obj = get_eleband_analysis(DIR_ELEBAND, DIR_ELEBAND_RESULTS, poscar_relaxed, extract_raw_data=True, extract_plot=True)
        elif i == PH:
            print('Now constructing phonon displacements and submitting jobs for force calculations.')
            ph_dir = checkPath(outDirName + PHONOPY_DIR_NAME)
            if not os.path.isdir(ph_dir):
                print(f"Making directory {ph_dir}...")
                mkdir(PHONOPY_DIR_NAME, outDirName)
            poscar_relaxed.write_file(ph_dir + POSCAR_UNIT_NAME)
            gsz = str(user_input_settings.get_super_dim())
            superdim_str = ' '.join((gsz, gsz, "1"))
            ph_preprocess(ph_dir, None, Poscar_unitcell_name=POSCAR_UNIT_NAME, 
                          onejob=False, user_input_settings=user_input_settings, 
                          supercellDim=superdim_str)
            print('Phonon displacement calculations kicked off successfully')
        else: 
            # Only case left is that we have phonon calculations to do.
            # First no matter what we need to preprocess to get FORCE_SETS. 
            if not ph_has_preprocessed:
                # Returns the directory name that all phonopy calculations are stored in.
                print('Preprocessing phonon analyses with phonopy displacements and VASP calculations...')
                incar_ph = copy.deepcopy(incar_selfcon)
                incar_ph = modifyIncar(incar_ph, addArr=[('SIGMA', SIGMA['narrow'])])
                DIR_PHONOPY = ph_prepare_for_analysis(outDirName, incar_ph, kpoints_mesh_ph, poscar_relaxed, potcar)
                DIR_PHONOPY = checkPath(DIR_PHONOPY)
                os.chdir(DIR_PHONOPY) # Go to that working directory
                print('Phonopy calculations subdirectory sent to postprocessor: %s'%(DIR_PHONOPY))
                ph_has_preprocessed = True
            else:
                print('Phonopy preprocessing already done. Proceeding directly to calculations...')
                os.chdir(DIR_PHONOPY) # Go to that working directory

            # Split into two analyses depending on whether band or dos
            if i == PHDOS:
                print('Now running phononic DOS calculations.')
                mkdir(PHDOS, DIR_PHONOPY)
                DIR_PHDOS = checkPath(DIR_PHONOPY + PHDOS)
                print('Conducting phonon total DOS analyses...')
                # Call the DOS analysis with the relevant parameters
                # NOTE: the poscar here is just for the atom names in the PUC so we don't need to do anything.
                ph_get_dos(kpoints_mesh_ph, poscar_relaxed, DIR_PHDOS, DIR_PHONOPY + POSCAR_UNIT_NAME)
            
            elif i == PHBAND:
                print('Now running phononic band structure calculations.')
                mkdir(PHBAND, DIR_PHONOPY)
                DIR_PHBAND = checkPath(DIR_PHONOPY + PHBAND)
                print('Conducting phonon band structure analyses...')
                # NOTE: the poscar here is just for the atom names in the PUC so we don't need to do anything.
                ph_get_band(kpoints_line, poscar_relaxed, DIR_PHBAND, DIR_PHONOPY + POSCAR_UNIT_NAME)

    
    # When it's all said and done, clean all the files written to the starter batch file directory
    cleanRelevantFiles(THIS_DIR)

    