from __directory_searchers import find # Will let us look for the different POSCAR-XYZ displacement files
from __dirModifications import move, copy, rm, mkdir

from ___constants_phonopy import *
from ___constants_names import *

import subprocess
import sys
from pymatgen.io.vasp.inputs import Poscar, VaspInput

# Function to handle all preprocessing of phonopy before the force calulation.
def phPreProcess(incar, kpoints_mesh, potcar, supercellDim=SUPER_DIM, Poscar_unitcell_name=POSCAR_UNIT_NAME, phDir=DIR_PHONOPY):
    # We need the h
    print(PHONOPY_DISP_MSG)

    try:
        CMD_GET_DISPLACEMENTS = ['phonopy', '-d', '--dim={}'.format(supercellDim), '-c', phDir + POSCAR_UNIT_NAME] # -c flag to indicate to phonopy where unit cell POSCAR is
        phStarterMsg = subprocess.run(CMD_GET_DISPLACEMENTS, capture_output=True, universal_newlines=True)
        print(phStarterMsg.stdout)

        # Check for the number of POSCAR-XYZ's made and organize them accordingly.
        poscarArray = find('POSCAR-*', DIR_PHONOPY)
        numPoscars = len(poscarArray)
        print('{} displacement files found.'.format(len(poscarArray)))
    except Exception as err:
        sys.exit('Phonopy preprocessing error: ' + err)

    if numPoscars == 0:
        sys.exit(PHONOPY_DISP_ERR_1)
    elif numPoscars > 999: # Obviously this should never happen but I'm all about caution with error handling
        sys.exit(PHONOPY_DISP_ERR_2)
    else:
        dispNum = 0 # Initialize outside loop block so we can return it
        for i in range(0, numPoscars):
            # Make the subfolders for each displacement and move the right POSCAR-XYZ there.
            dispNum = (poscarArray[i])[-3:] # Gives the XYZ of POSCAR-XYZ. THIS IS A STRING!
            try:
                newSubdir = PHDISP_DIR_NAME%(dispNum)
                makeNewSubdir = mkdir(newSubdir, DIR_PHONOPY)
                # NOTE: uncomment the next line and adjust the poscar input of the vasp object below if you want POSCAR-XYZ out of phonopy root dir
                # moveNewSubdir = move(poscarArray[i], DIR_PHONOPY, DIR_PHONOPY + newSubdir, POSCAR_NAME) # Rename it as well as move so we can just immediately run vasp

                # We also need to write the VASP input files into each of the new disp-XYZ folders so we can perform the VASP calculation
                poscar = Poscar.from_file(DIR_PHONOPY + poscarArray[i])
                vaspObj = VaspInput(incar, kpoints_mesh, poscar, potcar)
                # TODO: remove this write thing when VASP pymatgen gets working.
                vaspObj.write_input(output_dir=DIR_PHONOPY + newSubdir)
                
            except Exception as err:
                print(makeNewSubdir)
                # print(moveNewSubdir)
                sys.exit('Error in preprocessing phonopy (generating displacements and moving them to the right places): ' + err)

    return dispNum # Returns last displacement number so we can get force sets.

# Takes as input dispNum, the output of phPreProcess
def phGenerateForceSets(dispNum, dirName=DIR_PHONOPY):
    # We run a generator of force sets with phonopy on the shell
    phForce_output = subprocess.run('phonopy -f %s{001..%s}/%s'%(PHDISP_STATIC_NAME, dispNum, VASP_RUN_XML_NAME), shell=True, capture_output=True, universal_newlines=True)
    print(phForce_output.stdout)
    # TODO: If force sets created then proceed, if not exit with error. Use find function to find out.