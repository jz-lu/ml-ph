import sys
import subprocess

# Declare constants and hard strings
NO_INPUT_ERR_MSG = 'Usage: python3 _postProcess_relaxation.py <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband).'
BAD_INPUT_ERR_MSG = 'Error: invalid command line arguments.'
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND]
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full' # File path to the scripts

# Functions for editing the input files depending on the calculation we need.
def parseIncar():
    print('hi')

def parseKpoints():
    print('hi')

def organizePhonopyPoscars():
    print('hi')

def preprocessVaspForceCalc(type, )

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

    print(PHBAND)

# Process the command-line arguments
cmdargs = sys.argv
if len(cmdargs) == 1:
    sys.exit(NO_INPUT_ERR_MSG)
else:
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

            
    