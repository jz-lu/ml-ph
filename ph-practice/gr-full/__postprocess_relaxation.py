import sys, os
import subprocess
import matplotlib.pyplot as plt
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure

from __input_modifiers import modifyIncar, newKpoints, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy modules to the command line pipeline for basic linux commands
from __directory_searchers import checkPath
from __ph_processing import ph_prepare_for_analysis
from __run_vasp import run_vasp

from __get_eledos_analysis import get_eledos_analysis
from __get_eleband_analysis import get_eleband_analysis
from __get_elecombined_analysis import get_elecombined_analysis

from ___constants_vasp import NEDOS, ICHARG
from ___constants_names import *
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG


# Process the command-line arguments
def postProcess_relaxation(dirName, unrelaxed_vaspObj, calculation_list, kpoints_line): # vasp input object generated at beginning of relaxation
    dirName = checkPath(dirName)

    # Get the relevant objects for the next set of calculations
    chgcar = Chgcar.from_file(dirName + CHGCAR_NAME) # Need to retrieve the last charge density to get good calculations
    poscar_relaxed = Poscar.from_file(dirName + CONTCAR_NAME) # Get the relaxed positions from CONTCAR
    relaxation_incar = unrelaxed_vaspObj['INCAR'] # For the inspective code reviewer: these are pymatgen's hard strings so no need to collect them
    potcar = unrelaxed_vaspObj['POTCAR']
    
    # Incar changes for various calculations
    incar_selfcon = getSelfConNoRelIncar(relaxation_incar)
    incar_nonselfcon = getNonSelfConNoRelIncar(relaxation_incar)
    
    # Kpoints needs to be rebuilt since we want denser sampling for DOS calculations
    kpoints_mesh_nonrelax = newKpoints(dirName, 'mesh', poscar_relaxed)

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
            mkdir(i, dirName) # Create a subfolder for the analysis
            DIR_ELEDOS = checkPath(dirName + ELEDOS)
            mkdir(OUTPUT_DIR_NAME, DIR_ELEDOS) # subsubfolder for result storage
            DIR_ELEDOS_RESULTS = checkPath(DIR_ELEDOS + OUTPUT_DIR_NAME)

            # In addition to self-consistent calculation, we need to also add on the NEDOS parameter to sample the space for eledos
            incar_eledos = modifyIncar(incar_selfcon, addArr=[('NEDOS', NEDOS)])

            # We will need the standard VASP IO Object plus CHGCAR.
            eledos_vasp_obj = VaspInput(incar_eledos, kpoints_mesh_nonrelax, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})

            # Run vasp nonrelaxation, self-consistent
            run_vasp(eledos_vasp_obj, DIR_ELEDOS)

            # Get the analysis plots and data
            eledos_obj = get_eledos_analysis(DIR_ELEDOS, DIR_ELEDOS_RESULTS, poscar_relaxed, extract_raw_data=True, extract_plot=True)
            if os.path.isfile(DIR_ELEDOS + DOSCAR_NAME):
                move(DOSCAR_NAME, DIR_ELEDOS, DIR_ELEDOS_RESULTS)

        elif i == ELEBAND:
            mkdir(i, dirName) # Create a subfolder for the analysis
            DIR_ELEBAND = checkPath(dirName + ELEBAND)
            mkdir(OUTPUT_DIR_NAME, DIR_ELEBAND) # subsubfolder for result storage
            DIR_ELEBAND_RESULTS = checkPath(DIR_ELEBAND + OUTPUT_DIR_NAME)

            # We use the line kpoints file that we imported in the command line parsing start.py
            eleband_vasp_obj = VaspInput(incar_nonselfcon, kpoints_line, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})

            # Run vasp
            run_vasp(eleband_vasp_obj, DIR_ELEBAND)

            # Get the analysis plots and data
            eleband_obj = get_eleband_analysis(DIR_ELEBAND, DIR_ELEBAND_RESULTS, poscar_relaxed, extract_raw_data=True, extract_plot=True)

        else: 
            # Only case left is that we have phonon calculations to do
            # First no matter what we need to preprocess to get FORCE_SETS. 
            if not ph_has_preprocessed:
                # Returns the directory name that all phonopy calculations are stored in.
                DIR_PHONOPY = ph_prepare_for_analysis(dirName, incar_selfcon, kpoints_mesh_nonrelax, poscar_relaxed, potcar)
                DIR_PHONOPY = checkPath(DIR_PHONOPY)
                ph_has_preprocessed = True

            # Split into two analyses depending on whether band or dos
            if i == PHDOS:
                DIR_PHDOS = checkPath(DIR_PHONOPY + PHDOS)
            
            elif i == PHBAND:
                DIR_PHBAND = checkPath(DIR_PHONOPY + PHBAND)

            
    # Parse the full plot if there is one
    if (eleband_obj != None) and (eledos_obj != None):
        mkdir(COMBINED_ELE_OUTPUTS_NAME, dirName)
        combined_plot = get_elecombined_analysis(dirName + COMBINED_ELE_OUTPUTS_NAME, eledos_obj, eleband_obj)

    