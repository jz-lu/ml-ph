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

# Declare global checker variables
# hasBeenCalled = { 'getPhDisp': False,  }



# Organizes the generated displacement POSCAR-XYZ's from phonopy into nice subdirectories.
def organizePhonopyVaspOuts():
    # Algo outline: in a loop through XYZ make a copy of each input file into a disp-XYZ subfolder, and start a VASP job for  each
    # Define another function called cleanupPhonopyVaspOuts to get rid of extraneous inputs after the job is done, leaving just vasprun.xml.
    # Define another function to create FORCE_SETS up a folder (not in disp-XYZ subfolders). That's all before plotting!
    # In fact all of this should be wrapped up in a small function that does all the preprocessing before knowing whether
    # we want ph DOS or band, then just in the actual function do the plotting.
    # Make a checker bool global var to see if the processing was called already for a prev calc. If it was, then don't call it again.
    print('foo')
# Functions for processing. For each we first make a subfolder from where we are, in the original folder
# We then process the specific needs of the calculation

def processEledos():

    # For eledos, we need to copy the relaxed files, CONTCAR -> POSCAR,
    # modify INCAR to stop relaxation, use CHGCAR, and keep relaxation, 


    print(ELEDOS)

def processEleband():

    print(ELEBAND)

def processPhdos():

    print(PHDOS)

def processPhband():

    # TODO: move the proper input files in here for the processed vasp calculations, then run genDisp()

    print(PHBAND)

# Process the command-line arguments
cmdargs = sys.argv
if len(cmdargs) == 1:
    sys.exit(GENERAL_ERR_USAGE_MSG)
else:
    cmdargs = list(dict.fromkeys(cmdargs)) # converts to dict, which removes duplicates, and back to list
    del cmdargs[0]
    for i in cmdargs:
        if i not in CMD_LINE_ARG_LIST:
            print(GENERAL_ERR_USAGE_MSG)
            sys.exit(BAD_INPUT_ERR_MSG)
        # Makes sure that we do not have invalid command line arguments before we begin postprocessing

    # One more check to do: make sure that all the necessary VASP I/O files are there from the relaxation.
    # TODO: if vasp pmg runner is set up properly, make sure to change these getter directories to the output dir of the vasprun class
    original_incar = Incar.from_file(ROOT + INCAR_RELAXATION_NAME)
    incar_selfcon = getSelfConNoRelIncar(original_incar)
    incar_nonselfcon = getNonSelfConNoRelIncar(original_incar)
    
    # Besides INCAR no object needs changes between different analysis types
    poscar_relaxed = Poscar.from_file(ROOT + CONTCAR_NAME)
    potcar = Potcar.from_file(ROOT + POTCAR_NAME)
    chgcar = Chgcar.from_file(ROOT + CHGCAR_NAME)
    
    # Kpoints object in general may need an increase in density sampling depending on relaxation parameters
    original_kpoints_mesh = Kpoints.from_file(ROOT + KPOINTS_MESH_NAME)
    kpoints_mesh_relaxed = modifyKpoints('mesh')
    # Momentarily we wrap all these into a neat VaspInput class's object.

    # Parse command line and direct necessary function calls
    for i in cmdargs:
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

            
    