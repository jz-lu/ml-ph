import sys
import subprocess
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure

from __input_modifiers import modifyIncar, newKpoints, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy modules to the command line pipeline for basic linux commands
from __directory_searchers import checkPath
from __ph_processing import ph_preprocess, ph_generate_forcesets
from __run_vasp import run_vasp

from ___constants_vasp import NEDOS, ICHARG
from ___constants_names import *
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG


# Process the command-line arguments
def postProcess_relaxation(dirName, relaxed_vaspObj, calculation_list, kpoints_line): # vasp input object generated at beginning of relaxation
    dirName = checkPath(dirName)

    # Get the relevant objects for the next set of calculations
    chgcar = Chgcar.from_file(dirName + CHGCAR_NAME) # Need to retrieve the last charge density to get good calculations
    poscar_relaxed = Poscar.from_file(dirName + CONTCAR_NAME) # Get the relaxed positions from CONTCAR
    relaxation_incar = relaxed_vaspObj['INCAR'] # For the inspective code reviewer: these are pymatgen's hard strings so no need to collect them
    potcar = relaxed_vaspObj['POTCAR']
    
    # Incar changes for various calculations
    incar_selfcon = getSelfConNoRelIncar(relaxation_incar)
    incar_nonselfcon = getNonSelfConNoRelIncar(relaxation_incar)
    
    # Kpoints needs to be rebuilt since we want denser sampling for DOS calculations
    kpoints_mesh_nonrelax = newKpoints(dirName, 'mesh', poscar_relaxed)

    # Parse command line and direct necessary function calls
    for i in calculation_list:
        if i == ELEDOS:
            mkdir(i, dirName) # Create a subfolder for the analysis
            DIR_ELEDOS = dirName + ELEDOS + '/'

            # In addition to self-consistent calculation, we need to also add on the NEDOS parameter to sample the space for eledos
            incar_eledos = modifyIncar(incar_selfcon, addArr=[('NEDOS', NEDOS)])

            # We will need the standard VASP IO Object plus CHGCAR.
            eledos_vasp_obj = VaspInput(incar_eledos, kpoints_mesh_nonrelax, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})

            run_vasp(eledos_vasp_obj, DIR_ELEDOS)

        elif i == ELEBAND:
            mkdir(i, dirName) # Create a subfolder for the analysis
            DIR_ELEBAND = dirName + ELEBAND + '/'

            # We use the line kpoints file that we imported in the command line parsing start.py
            eleband_vasp_obj = VaspInput(incar_nonselfcon, kpoints_line, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})

            run_vasp(eleband_vasp_obj, DIR_ELEBAND)
        else: 
            # i.e. we have phonon calculations to do
            # First no matter what we need to preprocess to get FORCE_SETS. 
            # Only difference between band and DOS is the .conf file we generate.
            mkdir(PHONOPY_DIR_NAME, dirName)
            DIR_PHONOPY = dirName + PHONOPY_DIR_NAME + '/'

            # Due to the displacement invalidating charge densities, it is important that we use default charge densities to start
            # i.e. ICHARG = default = 2
            incar_selfcon_initChg = modifyIncar(incar_selfcon, addArr=[('ICHARG', ICHARG['default'])])

            # Note that in this object the poscar could be any valid poscar; we'll replace it in preprocessing by displacement poscar
            ph_preprocess_vasp_obj = VaspInput(incar_selfcon_initChg, kpoints_mesh_nonrelax, poscar_relaxed, potcar)
        
            # Run preprocessing
            lastDispNum = ph_preprocess(DIR_PHONOPY, ph_preprocess_vasp_obj) # returns a string with largest XYZ in POSCAR-XYZ for use in force sets generator
            
            # Generate force sets file
            ph_generate_forcesets(DIR_PHONOPY, lastDispNum)

            # TODO: conduct analysis on separate subpipes
            DIR_PHDOS = DIR_PHONOPY + PHDOS + '/'
            DIR_PHBAND = DIR_PHONOPY + PHBAND + '/'
            
    