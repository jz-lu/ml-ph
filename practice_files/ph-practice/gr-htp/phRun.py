import subprocess
import sys

# Declare constants and hard strings
RELATIVE_PH_DIR = '.'
POSCAR_UNIT = RELATIVE_PH_DIR + '/POSCAR_unit' # Unit cell POSCAR
DIM = '3 3 1'
BATCH_FILE = RELATIVE_PH_DIR + '/bat'
VASP_RUN_XML_NAME = 'vasprun.xml'

# Declare subprocess run arrays
CMD_GET_DISPLACEMENTS = ['phonopy', '-d', '--dim={}'.format(DIM), '-c', POSCAR_UNIT]

# A find function we'll need later for finding how many POSCARs there are
import os, fnmatch
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(name)
    return result

# First we run some things on terminal and compute the relaxation
# This is already done before as phonopy is  never called before preprocessed relaxation

# Get the displacements and organize them accordingly
# Input does not need to be checked as long as ele-DOS has exited successfully as they use same params
phStarterMsg = subprocess.run(CMD_GET_DISPLACEMENTS, capture_output=True, universal_newlines=True)
# Note that python automatically waits until the process is complete before moving to the next line of execution
# so no need to check. Returns a CompletedProcess object with the args if needed.
print('Starting phonopy calculations with displacement generation...\n{}'.format(phStarterMsg.stdout))

# Find all the POSCARs and perform a calculation for each one, then put vasp results into subdirectories
poscarArray = find('POSCAR-*', RELATIVE_PH_DIR)
numPoscars = len(poscarArray)

if numPoscars > 999:
    sys.exit("Too many POSCAR-displacement files to handle.") # This literally should never happen

pre = '' # Preface for POSCAR name, is it 00 or 0 or empty string?
for i in range(1, numPoscars + 1):
    if i <= 9:
        pre = '00'
    elif (10 <= i <= 99):
        pre = '0'
    else:
        pre = ''

    subprocess.run(['mkdir', '%s/disp-%s'%(RELATIVE_PH_DIR, pre+str(i))])

    # Run vasp with the specific displacement POSCAR and move that to the right subfolder
    print('Generated total of {} displacements'.format(numPoscars))
    print('Starting Vasp force calculations on phonopy generated displacements...')
    subprocess.run(['cp', '%s/POSCAR-%s'%(RELATIVE_PH_DIR, pre+str(i)), '%s/POSCAR'%(RELATIVE_PH_DIR)])
    phVaspJob = subprocess.run(['sbatch', BATCH_FILE], capture_output=True, universal_newlines=True)
    print(phVaspJob.stdout)

    curList = subprocess.run(['ls'], capture_output=True, universal_newlines=True)
    print('Directory list:', curList.stdout)

    subfolder = RELATIVE_PH_DIR + ('/disp-%s'%(pre+str(i)))
    subprocess.run(['mv', RELATIVE_PH_DIR+'/'+VASP_RUN_XML_NAME, subfolder])
    print('Job complete for displacement %s, run file stored in %s'%(pre+str(i), subfolder))

lastDisplacement = pre + str(numPoscars)
print('Last displacement:', lastDisplacement)
forceObj = subprocess.run(['phonopy', '-f', 'disp-{001..%s}'%(lastDisplacement), '-c', POSCAR_UNIT], capture_output=True, universal_newlines=True)
print(forceObj.stdout)
