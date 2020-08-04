from __directory_searchers import find # Will let us look for the different POSCAR-XYZ displacement files
from __dirModifications import move, copy, rm, mkdir

from ___constants_phonopy import *
from ___constants_names import *

import subprocess
import sys

def phPreProcess(supercellDim=SUPER_DIM, Poscar_unitcell_name=POSCAR_UNIT, phDir=PHONOPY_DIR_NAME):
    print(PHONOPY_DISP_MSG)
    CMD_GET_DISPLACEMENTS = ['phonopy', '-d', '--dim={}'.format(supercellDim), '-c', PH_ROOT + '/' + POSCAR_UNIT]
    phStarterMsg = subprocess.run(CMD_GET_DISPLACEMENTS, capture_output=True, universal_newlines=True)
    print(phStarterMsg.stdout)

    # Check for the number of POSCAR-XYZ's made and organize them accordingly.
    poscarArray = find('POSCAR-*', phDir)
    numPoscars = len(poscarArray)
    print('{} displacement files found.'.format(len(poscarArray)))

    if numPoscars == 0:
        sys.exit(PHONOPY_DISP_ERR_1)
    elif numPoscars > 999: # Obviously this should never happen but I'm all about caution with error handling
        sys.exit(PHONOPY_DISP_ERR_2)
    else:
        for i in range(0, numPoscars):
            # Make the subfolders for each displacement and move the right POSCAR-XYZ there.
            dispNum = (poscarArray[i])[-3:]
            newSubdir = 'disp-%s'%(dispNum)
            makeNewSubdir = mkdir(newSubdir, ROOT)
            moveNewSubdir = move(poscarArray[i], ROOT, ROOT + '/' + newSubdir, 'POSCAR') # Rename it so we can just immediately run vasp
            # subprocess.run(['mv', ROOT + '/' + poscarArray[i], ROOT + '/' + newSubdir])

            if (makeNewSubdir != '') or (moveNewSubdir != ''):
                sys.exit('Error in moving folders. Logs: %s\n%s'%(makeNewSubdir, moveNewSubdir))

            # Move all the other relevant files as well for VASP force calculations. Then in the post process of phonopy do a cleanup.
            
    return