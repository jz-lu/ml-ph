import sys
import subprocess
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure

from __input_modifiers import modifyIncar, newKpoints, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy modules to the command line pipeline for basic linux commands
from __directory_searchers import checkPath
from __ph_processing import ph_preprocess
from __run_vasp import run_vasp

from ___constants_vasp import NEDOS
from ___constants_names import *
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG


# Process the command-line arguments
def postProcess_relaxation(relaxed_vaspObj, calculation_list, dirName): # vasp input object generated at beginning of relaxation
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

    # Declare necessary directory strings. Note that exceot for the batch file path, all must end in '/' for postprocessing purposes
    DIR_ELEBAND = dirName + ELEBAND + '/'
    DIR_PHONOPY = dirName + PHONOPY_DIR_NAME + '/'
    DIR_PHDOS = DIR_PHONOPY + PHDOS + '/'
    DIR_PHBAND = DIR_PHONOPY + PHBAND + '/'

    # Parse command line and direct necessary function calls
    for i in calculation_list:
        if i == ELEDOS:
            mkdir(i, dirName) # Create a subfolder for the analysis
            DIR_ELEDOS = dirName + ELEDOS + '/'

            # In addition to self-consistent calculation, we need to also add on the NEDOS parameter to sample the space for eledos
            incar_eledos = modifyIncar(incar_selfcon, addArr=[('NEDOS', NEDOS)])

            # We will need the standard VASP IO Object plus CHGCAR.
            eledos_vasp_obj = VaspInput(incar_eledos, kpoints_mesh_nonrelax, poscar_relaxed, potcar, {CHGCAR_NAME: chgcar})

            # Now we run vasp
            run

            # TODO: change the copy files thing to a more robust VaspInput.from_directory. SAME FOR BAND.
        elif i == ELEBAND:
            mkdir(i, dirName) # Create a subfolder for the analysis

            # TODO: same TODO as ELEDOS but for band. CHGCAR into Vasp input object.

            try:
                # TODO: change the line to autogeneration as well.
                kpoints_line = Kpoints.from_file(dirName + KPOINTS_LINE_NAME)
            except Exception as err:
                print('Error:', err)
                sys.exit('Suggested source of error: you asked for a band structure calculation but did not give a (valid) line kpoints file named {}'.format(KPOINTS_LINE_NAME))

            # We will need the standard VASP IO Object plus CHGCAR.
            stdVaspObj = VaspInput(incar_nonselfcon, kpoints_line, poscar_relaxed, potcar)
            copy(CHGCAR_NAME, dirName, DIR_ELEBAND)
            stdVaspObj.write_input(DIR_ELEDOS)
        else: # i.e. we have phonon calculations to do
            # First no matter what we need to preprocess to get FORCE_SETS. 
            # Only difference between band and DOS is the .conf file we generate.
            mkdir(PHONOPY_DIR_NAME, dirName)
            stdVaspObj = VaspInput(incar_selfcon, kpoints_mesh_relaxed, poscar_relaxed, potcar)
            # TODO: is this necessary? Should phonopy run with ICHARG = 1 or default?
            copy(CHGCAR_NAME, dirName, DIR_PHONOPY)
            # TODO: if this is necessary, then figure out a way to put it into the vasp object container or pass it into phPreprocess()
            

            # TODO: phPreprocess()

            
    