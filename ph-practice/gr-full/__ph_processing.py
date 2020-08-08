from __directory_searchers import find # Will let us look for the different POSCAR-XYZ displacement files
from __dirModifications import move, copy, rm, mkdir
from __directory_searchers import checkPath, findFilesInDir, filesInDir
from __run_vasp import run_vasp
from __input_modifiers import modifyIncar
from __ph_movers import moveRelevantFiles

from ___constants_phonopy import *
from ___constants_names import *
from ___constants_misc import ERR_PH_FORCE_SETS, ERR_PH_FORCE_SETS_NOT_FOUND
from ___constants_vasp import *

import subprocess
import sys, os
from pymatgen.io.vasp.inputs import Poscar, VaspInput

# Function to handle all preprocessing of phonopy before the force calulation.
def ph_preprocess(dirName, vaspObj, supercellDim=SUPER_DIM, Poscar_unitcell_name=POSCAR_UNIT_NAME):
    print(PHONOPY_DISP_MSG)
    dirName = checkPath(dirName)

    try:
        CMD_GET_DISPLACEMENTS = ['phonopy', '-d', '--dim={}'.format(supercellDim), '-c', dirName + Poscar_unitcell_name] # -c flag to indicate to phonopy where unit cell POSCAR is
        phStarterMsg = subprocess.run(CMD_GET_DISPLACEMENTS, capture_output=True, universal_newlines=True)
        print(phStarterMsg.stdout)

        # Phonopy places all their files in the directory of this script. 
        # We want it in phDir, so we scane for all files that are not .py scripts or related and move them over.
        moveRelevantFiles(dirName)

        # Check for the number of POSCAR-XYZ's made and organize them accordingly.
        poscarArray = findFilesInDir(dirName, 'POSCAR-', 'start') # Again, phonopy hard string 'POSCAR-XYZ', no need to collect
        poscarArray.sort() # Make sure displacements in order
        numPoscars = len(poscarArray)
        print('Result: {} displacement files found.'.format(len(poscarArray)))
    except Exception as err:
        sys.exit('Phonopy preprocessing error: ' + err)

    if numPoscars == 0:
        sys.exit(PHONOPY_DISP_ERR_1)
    elif numPoscars > 999: # Obviously this should never happen but I'm all about caution with error handling
        sys.exit(PHONOPY_DISP_ERR_2)
    else:
        dispNum = '' # Initialize outside loop block so we can return it
        for i in range(0, numPoscars):
            # Make the subfolders for each displacement and move the right POSCAR-XYZ there.
            dispNum = (poscarArray[i])[-3:] # Gives the XYZ of POSCAR-XYZ. THIS IS A STRING!
            try:
                newSubdirName = PHDISP_DIR_NAME%(dispNum)
                mkdir(newSubdirName, dirName)

                # Perform the VASP calculation on this 
                # NOTE: no chgcar, the charge densities would be inaccurate now that we have displacements
                displaced_poscar = Poscar.from_file(dirName + poscarArray[i])
                ph_vasp_obj = VaspInput(vaspObj['INCAR'], vaspObj['KPOINTS'], displaced_poscar, vaspObj['POTCAR'])
                print('Starting VASP nonrelaxation forces calculation for displacement {}...'.format(dispNum))
                run_vasp(ph_vasp_obj, dirName + newSubdirName)
                print('VASP calculation complete.')
                
            except Exception as err:
                print('An error occurred while processing displacement POSCAR file {}'.format(dirName + poscarArray[i]))
                sys.exit('Error in preprocessing phonopy (parsing displacement files and running VASP force calculations): ' + err)

        print('Total number of displacement files generated: ' + dispNum)

    return dispNum # Returns last displacement number so we can get force sets.

# Takes as input dispNum, the output of phPreProcess
def ph_generate_forcesets(dirName, dispNum):
    dirName = checkPath(dirName)
    try:
        if dispNum ==  '001':
            print('Generating force sets now...')
            print('Running command to shell: phonopy -f %s001/%s'%(PHDISP_STATIC_NAME, VASP_RUN_XML_NAME))
            phForce_output = subprocess.run('phonopy -f %s001/%s'%(PHDISP_STATIC_NAME, VASP_RUN_XML_NAME), shell=True, capture_output=True, universal_newlines=True)
            print(phForce_output.stdout)
        else:
            # We run a generator of force sets with phonopy on the shell
            print('Generating force sets now...')
            print('Running command to shell: phonopy -f %s{001..%s}/%s'%(PHDISP_STATIC_NAME, dispNum, VASP_RUN_XML_NAME))
            phForce_output = subprocess.run('phonopy -f %s{001..%s}/%s'%(PHDISP_STATIC_NAME, dispNum, VASP_RUN_XML_NAME), shell=True, capture_output=True, universal_newlines=True)
            print(phForce_output.stdout)
    except Exception as err:
        print(ERR_PH_FORCE_SETS)
        print('Command log: {}'.format(phForce_output.stderr))
        sys.exit('Error:', err)

    # If force sets created then proceed, if not exit with error. Use find function to find out.
    if not os.path.isfile(THIS_DIR + PH_FORCE_SETS_NAME):
        sys.exit(ERR_PH_FORCE_SETS_NOT_FOUND)
    
    # Move FORCE_SETS to the proper place
    move(PH_FORCE_SETS_NAME, THIS_DIR, dirName)
    print(PH_FORCE_SETS_NAME + ' successfully moved to ' + dirName + '. Phonopy postprocessing complete.')

    return 
    
# Get a neat function that uses the above functions to handle all processing
def ph_prepare_for_analysis(rootDirName, incar_selfcon, kpoints_mesh_nonrelax, poscar_relaxed, potcar):
    rootDirName = checkPath(rootDirName)
    mkdir(PHONOPY_DIR_NAME, rootDirName)
    DIR_PHONOPY = checkPath(rootDirName + PHONOPY_DIR_NAME)

    # Write out the poscar into the phonopy directory we made, so we can use it for calling phonopy
    try:
        poscar_relaxed.write_file(DIR_PHONOPY + POSCAR_UNIT_NAME)
    except Exception as err:
        print('Error in writing the %s file to %s.'%(POSCAR_UNIT_NAME, DIR_PHONOPY))
        sys.exit(err)

    # Due to the displacement invalidating charge densities, it is important that we use default charge densities to start
    # i.e. ICHARG = default = 2
    incar_selfcon_initChg = modifyIncar(incar_selfcon, addArr=[('ICHARG', ICHARG['default'])])

    # Note that in this object the poscar could be any valid poscar; we'll replace it in preprocessing by displacement poscar
    ph_preprocess_vasp_obj = VaspInput(incar_selfcon_initChg, kpoints_mesh_nonrelax, poscar_relaxed, potcar)

    # Run preprocessing
    lastDispNum = ph_preprocess(DIR_PHONOPY, ph_preprocess_vasp_obj) # returns a string with largest XYZ in POSCAR-XYZ for use in force sets generator
    
    # Generate force sets file
    ph_generate_forcesets(DIR_PHONOPY, lastDispNum)

    return DIR_PHONOPY

    