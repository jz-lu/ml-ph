import sys
import subprocess
from pymatgen.io.vasp.inputs import Incar, Kpoints
from pymatgen.core.structure import Structure

from __input_modifiers import modifyIncar, modifyKpoints # Input modifiers
from __dirModifications import move, copy, mkdir, rm # Easy "APIs" to the command line pipeline for basic linux commands

# Declare constants
NO_INPUT_ERR_MSG = 'Usage: python3 _postProcess_relaxation.py <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband).'
BAD_INPUT_ERR_MSG = 'Error: invalid command line arguments.'
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND]
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full' # File path to the scripts, must be in same folder as relaxation.

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

# Functions for processing. For each we first make a subfolder from where we are, in the original folder
# We then process the specific needs of the calculation

def processEledos():
    dirPath = ROOT + '/' + ELEDOS
    subprocess.run(['mkdir', dirPath])

    # For eledos, we need to copy the relaxed files, CONTCAR -> POSCAR,
    # modify INCAR to stop relaxation, 
    print(ELEDOS)

def processEleband():
    dirPath = ROOT + '/' + ELEBAND
    subprocess.run(['mkdir', dirPath])

    print(ELEBAND)

def processPhdos():
    dirPath = ROOT + '/' + PHDOS
    subprocess.run(['mkdir', dirPath])

    print(PHDOS)

def processPhband():
    dirPath = ROOT + '/' + PHBAND
    subprocess.run(['mkdir', dirPath])

    # TODO: move the proper input files in here for the processed vasp calculations, then run genDisp()

    print(PHBAND)

# Process the command-line arguments
cmdargs = sys.argv
if len(cmdargs) == 1:
    sys.exit(NO_INPUT_ERR_MSG)
else:
    subprocess.run(['cp', ROOT + '/CONTCAR', ROOT + '/' + POSCAR_UNIT]) # Get the relaxed unit cell POSCAR relaxed

    # Parse command line and direct necessary function calls
    del cmdargs[0]
    isCreated = {
            ELEDOS: False,
            ELEBAND: False,
            PHDOS: False,
            PHBAND: False
        } # Ensures we don't create the same thing over and over if user specifies same flag multiple times

    # Parse the command line arguments accordingly
    for i in cmdargs:
        if i not in CMD_LINE_ARG_LIST:
            print(BAD_INPUT_ERR_MSG)
            sys.exit(NO_INPUT_ERR_MSG)
        else:
            if not isCreated[i]:
                if i == ELEDOS:
                    processEledos()
                elif i == ELEBAND:
                    processEleband()
                elif i == PHDOS:
                    processPhdos()
                else:
                    processPhband()
                isCreated[i] = True

            
    