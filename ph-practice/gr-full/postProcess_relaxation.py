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
def postProcess_relaxation(vaspObj, calculation_list, vaspRelaxationOutDir=ROOT): # vasp input object generated at beginning of relaxation

    # Make sure that all the necessary VASP I/O files are there from the relaxation.
    # TODO: if vasp pmg runner is set up properly, make sure to change these getter directories to the output dir of the vasprun class
    original_incar = vaspObj['INCAR'] # These are pymatgen's hard strings so no need to collect them
    incar_selfcon = getSelfConNoRelIncar(original_incar)
    incar_nonselfcon = getNonSelfConNoRelIncar(original_incar)
    
    # Besides INCAR no object needs changes between different analysis types
    poscar_relaxed = vaspObj['POSCAR']
    potcar = vaspObj['POTCAR']
    chgcar = Chgcar.from_file(ROOT + CHGCAR_NAME) # Only one we need to retrieve
    
    # Kpoints object in general may need an increase in density sampling depending on relaxation parameters
    original_kpoints_mesh = Kpoints.from_file(ROOT + KPOINTS_MESH_NAME)
    kpoints_mesh_relaxed = modifyKpoints('mesh')
    # Momentarily we wrap all these into a neat VaspInput class's object.

    # Parse command line and direct necessary function calls
    for i in calculation_list:
        if i == ELEDOS:
            mkdir(i, ROOT) # Create a subfolder for the analysis

            # We will need the standard VASP IO Object plus CHGCAR.
            stdVaspObj = VaspInput(incar_selfcon, kpoints_mesh_relaxed, poscar_relaxed, potcar)

            # TODO: delete the next line and put CHGCAR as an object into the vasp obj when you figure out how.
            copy(CHGCAR_NAME, ROOT, DIR_ELEDOS)
            stdVaspObj.write_input(DIR_ELEDOS)

            # TODO: change the copy files thing to a more robust VaspInput.from_directory. SAME FOR BAND.
        elif i == ELEBAND:
            mkdir(i, ROOT) # Create a subfolder for the analysis

            # TODO: same TODO as ELEDOS but for band. CHGCAR into Vasp input object.

            try:
                # TODO: change the line to autogeneration as well.
                kpoints_line = Kpoints.from_file(ROOT + KPOINTS_LINE_NAME)
            except Exception as err:
                print('Error:', err)
                sys.exit('Suggested source of error: you asked for a band structure calculation but did not give a (valid) line kpoints file named {}'.format(KPOINTS_LINE_NAME))

            # We will need the standard VASP IO Object plus CHGCAR.
            stdVaspObj = VaspInput(incar_nonselfcon, kpoints_line, poscar_relaxed, potcar)
            copy(CHGCAR_NAME, ROOT, DIR_ELEBAND)
            stdVaspObj.write_input(DIR_ELEDOS)
        else: # i.e. we have phonon calculations to do
            # First no matter what we need to preprocess to get FORCE_SETS. 
            # Only difference between band and DOS is the .conf file we generate.
            mkdir(PHONOPY_DIR_NAME, ROOT)
            stdVaspObj = VaspInput(incar_selfcon, kpoints_mesh_relaxed, poscar_relaxed, potcar)
            # TODO: is this necessary? Should phonopy run with ICHARG = 1 or default?
            copy(CHGCAR_NAME, ROOT, DIR_PHONOPY)
            # TODO: if this is necessary, then figure out a way to put it into the vasp object container or pass it into phPreprocess()
            

            # TODO: phPreprocess()

            
    