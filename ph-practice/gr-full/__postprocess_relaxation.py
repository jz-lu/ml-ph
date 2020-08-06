import sys
import subprocess
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure

from __input_modifiers import modifyKpoints, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy "APIs" to the command line pipeline for basic linux commands
from ___constants_names import *
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG

from __phProcessing import phPreProcess

# Process the command-line arguments
def postProcess_relaxation(relaxed_vaspObj, calculation_list, dirName): # vasp input object generated at beginning of relaxation
    # Get the relevant objects for the next set of calculations
    chgcar = Chgcar.from_file(dirName + CHGCAR_NAME) # Need to retrieve the last charge density to get good calculations
    relaxed_poscar = Poscar.from_file(dirName + CONTCAR_NAME) # Get the relaxed positions
    relaxation_incar = relaxed_vaspObj['INCAR'] # These are pymatgen's hard strings so no need to collect them
    potcar = relaxed_vaspObj['POTCAR']
    
    # Incar changes for various calculations
    incar_selfcon = getSelfConNoRelIncar(relaxation_incar)
    incar_nonselfcon = getNonSelfConNoRelIncar(relaxation_incar)
    
    # Kpoints object in general may need an increase in density sampling depending on relaxation parameters
    original_kpoints_mesh = Kpoints.from_file(dirName + KPOINTS_MESH_NAME)
    kpoints_mesh_relaxed = modifyKpoints(dirName, 'mesh', relaxed_poscar)
    # Momentarily we wrap all these into a neat VaspInput class's object.

    # Declare necessary directory strings. Note that exceot for the batch file path, all must end in '/' for postprocessing purposes
    DIR_ELEDOS = dirName + ELEDOS + '/'
    DIR_ELEBAND = dirName + ELEBAND + '/'
    PHDISP_DIR_NAME = PHDISP_STATIC_NAME + '%s/'
    DIR_PHONOPY = dirName + PHONOPY_DIR_NAME + '/'
    DIR_PHDOS = DIR_PHONOPY + PHDOS + '/'
    DIR_PHBAND = DIR_PHONOPY + PHBAND + '/'

    # Parse command line and direct necessary function calls
    for i in calculation_list:
        if i == ELEDOS:
            mkdir(i, dirName) # Create a subfolder for the analysis

            # We will need the standard VASP IO Object plus CHGCAR.
            stdVaspObj = VaspInput(incar_selfcon, kpoints_mesh_relaxed, poscar_relaxed, potcar)

            # TODO: delete the next line and put CHGCAR as an object into the vasp obj when you figure out how.
            copy(CHGCAR_NAME, dirName, DIR_ELEDOS)
            stdVaspObj.write_input(DIR_ELEDOS)

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

            
    