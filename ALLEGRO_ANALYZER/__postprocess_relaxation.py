import os
import subprocess
import matplotlib.pyplot as plt
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure

from __input_modifiers import modifyIncar, modifyKpointsMesh, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy modules to the command line pipeline for basic linux commands
from __directory_searchers import checkPath
from __run_vasp import run_vasp
from __ph_processing import ph_prepare_for_analysis
from __get_ph_analysis import ph_get_dos, ph_get_band

from __get_eledos_analysis import get_eledos_analysis
from __get_eleband_analysis import get_eleband_analysis
from __get_elecombined_analysis import get_elecombined_analysis

from ___constants_vasp import NEDOS, ICHARG
from ___constants_names import *
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG
from ___constants_phonopy import POSCAR_UNIT_NAME


# Process the command-line arguments
def postProcess_relaxation(outDirName, relaxation_dirName, unrelaxed_vaspObj, calculation_list, kpoints_line): # vasp input object generated at beginning of relaxation
    print('Postprocessing VASP relaxation beginning...')
    outDirName = checkPath(outDirName)
    relaxation_dirName = checkPath(relaxation_dirName)

    # Get the relevant objects for the next set of calculations
    chgcar = Chgcar.from_file(relaxation_dirName + CHGCAR_NAME) # Need to retrieve the last charge density to get good calculations
    poscar_relaxed = Poscar.from_file(relaxation_dirName + CONTCAR_NAME) # Get the relaxed positions from CONTCAR
    relaxation_incar = unrelaxed_vaspObj['INCAR'] # For the inspective code reviewer: these are pymatgen's hard strings so no need to collect them
    potcar = unrelaxed_vaspObj['POTCAR']
    
    # Incar changes for various calculations
    incar_selfcon = getSelfConNoRelIncar(relaxation_incar)
    incar_nonselfcon = getNonSelfConNoRelIncar(relaxation_incar)
    
    # Kpoints needs to be modified to have a denser sampling for nonrelaxation calculations -> more precision
    kpoints_mesh_nonrelax = unrelaxed_vaspObj['KPOINTS']
    kpoints_mesh_nonrelax = modifyKpointsMesh(kpoints_mesh_nonrelax)

    # Objects useful for plotting data
    eledos_obj = None
    eleband_obj = None
    combinedPlot = None
    
    # Since all processing for phonopy is done at once, we should flag so that if we want phband and phdos we don't preprocess twice
    ph_has_preprocessed = False
    DIR_PHONOPY = None # This will be set during first processing

    # Parse command line and direct necessary function calls
    for i in calculation_list:
        if i == ELEDOS:
            print('Now running electronic DOS calculations.')
            mkdir(i, outDirName) # Create a subfolder for the analysis
            DIR_ELEDOS = checkPath(outDirName + ELEDOS)
            mkdir(OUTPUT_DIR_NAME, DIR_ELEDOS) # subsubfolder for result storage
            DIR_ELEDOS_RESULTS = checkPath(DIR_ELEDOS + OUTPUT_DIR_NAME)

            # In addition to self-consistent calculation, we need to also add on the NEDOS parameter to sample the space for eledos
            incar_eledos = modifyIncar(incar_selfcon, addArr=[('NEDOS', NEDOS)])

            # We will need the standard VASP IO Object plus CHGCAR.
            eledos_vasp_obj = VaspInput(incar_eledos, kpoints_mesh_nonrelax, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})

            print('Running VASP to get data for electronic DOS...')
            # Run vasp nonrelaxation, self-consistent
            run_vasp(eledos_vasp_obj, DIR_ELEDOS)

            # Get the analysis plots and data
            print('Parsing VASP run and retrieving DOS data...')
            eledos_obj = get_eledos_analysis(DIR_ELEDOS, DIR_ELEDOS_RESULTS, poscar_relaxed, extract_raw_data=True, extract_plot=True)
            if os.path.isfile(DIR_ELEDOS + DOSCAR_NAME):
                move(DOSCAR_NAME, DIR_ELEDOS, DIR_ELEDOS_RESULTS)

        elif i == ELEBAND:
            print('Now running electronic band structure calculations.')
            mkdir(i, outDirName) # Create a subfolder for the analysis
            DIR_ELEBAND = checkPath(outDirName + ELEBAND)
            mkdir(OUTPUT_DIR_NAME, DIR_ELEBAND) # subsubfolder for result storage
            DIR_ELEBAND_RESULTS = checkPath(DIR_ELEBAND + OUTPUT_DIR_NAME)

            # We use the line kpoints file that we imported in the command line parsing start.py
            eleband_vasp_obj = VaspInput(incar_nonselfcon, kpoints_line, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})
            
            print('Running VASP for electronic band calculations...')
            
            # Run vasp
            run_vasp(eleband_vasp_obj, DIR_ELEBAND)

            # Get the analysis plots and data
            print('Parsing VASP run and retrieving band structure data...')
            eleband_obj = get_eleband_analysis(DIR_ELEBAND, DIR_ELEBAND_RESULTS, poscar_relaxed, extract_raw_data=True, extract_plot=True)

        else: 
            # Only case left is that we have phonon calculations to do
            # First no matter what we need to preprocess to get FORCE_SETS. 
            if not ph_has_preprocessed:
                # Returns the directory name that all phonopy calculations are stored in.
                print('Preprocessing phonon analyses with phonopy displacements and VASP calculations...')
                DIR_PHONOPY = ph_prepare_for_analysis(outDirName, incar_selfcon, kpoints_mesh_nonrelax, poscar_relaxed, potcar)
                DIR_PHONOPY = checkPath(DIR_PHONOPY)
                ph_has_preprocessed = True

            # Split into two analyses depending on whether band or dos
            if i == PHDOS:
                print('Now running phononic DOS calculations.')
                mkdir(PHDOS, DIR_PHONOPY)
                DIR_PHDOS = checkPath(DIR_PHONOPY + PHDOS)
                print('Conducting phonon total DOS analyses...')
                # Call the DOS analysis with the relevant parameters
                # NOTE: the poscar here is just for the atom names in the PUC so we don't need to do anything.
                ph_get_dos(kpoints_mesh_nonrelax, poscar_relaxed, DIR_PHDOS, DIR_PHONOPY + POSCAR_UNIT_NAME)
            
            elif i == PHBAND:
                print('Now running phononic band structure calculations.')
                mkdir(PHBAND, DIR_PHONOPY)
                DIR_PHBAND = checkPath(DIR_PHONOPY + PHBAND)
                print('Conducting phonon band structure analyses...')
                # NOTE: the poscar here is just for the atom names in the PUC so we don't need to do anything.
                ph_get_band(kpoints_line, poscar_relaxed, DIR_PHBAND, DIR_PHONOPY + POSCAR_UNIT_NAME)

            
    # Parse the full plot if there is one
    if (eleband_obj != None) and (eledos_obj != None):
        print('Both electronic DOS and band structure calculated, combined plot available.')
        print('Preparing a combined electron DOS and band structure plot...')
        mkdir(COMBINED_ELE_OUTPUTS_NAME, outDirName)
        combined_plot = get_elecombined_analysis(outDirName + COMBINED_ELE_OUTPUTS_NAME, eledos_obj, eleband_obj)
    
    print(PRGM_END_CARD)

    