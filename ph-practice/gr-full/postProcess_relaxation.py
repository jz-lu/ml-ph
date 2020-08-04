import sys
import subprocess
from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar, Potcar, VaspInput
from pymatgen.io.vasp.outputs import Chgcar
from pymatgen.core.structure import Structure

from __input_modifiers import modifyKpoints, getSelfConNoRelIncar, getNonSelfConNoRelIncar # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy "APIs" to the command line pipeline for basic linux commands
from ___constants_names import *
from ___constants_misc import BAD_INPUT_ERR_MSG, GENERAL_ERR_USAGE_MSG

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
    mkdir(ELEBAND, ROOT)

    # For eledos, we need to copy the relaxed files, CONTCAR -> POSCAR,
    # modify INCAR to stop relaxation, use CHGCAR, and keep relaxation, 
        

    print(ELEDOS)

def processEleband():
    mkdir(ELEBAND, ROOT)

    print(ELEBAND)

def processPhdos():
    mkdir(PHDOS, ROOT)

    print(PHDOS)

def processPhband():
    mkdir(PHBAND, ROOT)

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
    relaxed_poscar = Poscar.from_file(ROOT + CONTCAR_NAME)
    potcar = Potcar.from_file(ROOT + POTCAR_NAME)
    chgcar = Chgcar.from_file(ROOT + CHGCAR_NAME)
    
    # Kpoints object in general may need an increase in density sampling depending on relaxation parameters
    original_kpoints_mesh = Kpoints.from_file(ROOT + KPOINTS_MESH_NAME)
    kpoints_mesh_relaxed = modifyKpoints('mesh')
    # Momentarily we wrap all these into a neat VaspInput class's object.

    # Parse command line and direct necessary function calls

    for i in cmdargs:
        if i == ELEDOS:
            # We will need the standard VASP IO Object plus CHGCAR.
            processEledos()
        elif i == ELEBAND:
            processEleband()
        elif i == PHDOS:
            processPhdos()
        elif i == PHBAND:
            processPhband()

            
    