from ___constants_misc import *
from ___constants_names import *
from ___constants_vasp import *
from __directory_searchers import find
from __dirModifications import move, cat, checkPath
from __query_inputs import getPoscarObj, getInputName
import sys

from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints, Potcar

# Returns pymatgen Kpoints object
def buildRelaxationKpoints(dirName=ROOT, grid=RELAXATION_GRID_DENSITY, shift=RELAXATION_GRID_SHIFT, isRelaxed=False, writeOut=True):
    kpointsExists = find(KPOINTS_MESH_NAME, ROOT) # Look for existing KPOINTS
    if kpointsExists:
        try:
            kpoints = Kpoints.from_file(checkPath(dirName) + KPOINTS_MESH_NAME)
            print('KPOINTS file inputted by user. Importing input for computation...')
        except Exception as err:
            print('Error:', err)
            print('Suggestion: your user-inputted KPOINTS file is likely invalid. Check it')
            sys.exit()
    else:
        if grid == RELAXATION_GRID_DENSITY and shift == RELAXATION_GRID_SHIFT:
            print('No KPOINTS file found. Generating default relaxation mesh...')
        else:
            print('No KPOINTS file found. Generating mesh according to user inputted grid...')
        kpoints = Kpoints.gamma_automatic(grid, shift) # Only supports Gamma centered Monkhorst Pack for now
    
    kpoints.comment = getInputName() # Standardize labeling style

    if (not kpointsExists) and writeOut:
        print('Writing new KPOINTS file...')
        kpoints.write_file(dirName + KPOINTS_MESH_NAME)
    
    if not writeOut:
        print('KPOINTS object created without file writing.')

    return kpoints

# Returns pymatgen Incar object
def buildRelaxationIncar(dirName=ROOT, vdW=False, writeOut=True):
    # First let's build an array of default keys, which is helpful later
    default_keys = list(INCAR_DEFAULT_SETTINGS.keys())
    default_values = list(INCAR_DEFAULT_SETTINGS.values())
    vdW_keys = list(INCAR_VDW_SETTINGS.keys())
    vdW_values = list(INCAR_VDW_SETTINGS.values())

    incarExists = find(INCAR_RELAXATION_NAME, ROOT) # Look for existing INCAR
    if not incarExists:
        # Create a default
        print('No INCAR input found. Generating default relaxation settings...')
        incar = Incar.from_string('')
        for i in range(0, len(default_keys)):
            incar[default_keys[i]] = default_values[i]
    else:
        print('INCAR input found. Using parameters given and adding any other necessary ones...')
        try:
            incar = Incar.from_file(ROOT + INCAR_RELAXATION_NAME)
        except Exception as error:
            print('Error:', error)
            print('Suggestion: you probably have a file labeled INCAR in the specified directory that is invalid.')

        for i in range(0, len(default_keys)):
            if default_keys[i] not in incar:
                print('{} not in INCAR. Adding it automatically for relaxation calculations...'.format(default_keys[i]))
                incar[default_keys[i]] = default_values[i]
    
    if vdW:
        print('Adding van der Waals interaction parameters to INCAR...')
        for i in range(0, len(vdW_keys)):
            incar[vdW_keys[i]] = vdW_values[i]
    
    # Get a SYSTEM tag for the name, as a comment
    incar['SYSTEM'] = getInputName()

    if writeOut:
        print('Writing INCAR to file...')
        incar.write_file(ROOT + INCAR_RELAXATION_NAME)
    
    return incar

# Writes a properly constructed POTCAR (not an object), use only if pymatgen method below isn't working!
# Returns 0 on success, 1 on failure
def buildPotcar_noPMG(dirName=ROOT, potcarDir=POT_NOPMG_DIR, poscarObj=getPoscarObj()):
    atoms = poscarObj.site_symbols
    # print(atoms)

    dirName = checkPath(dirName)
    potcarDir = checkPath(potcarDir)

    fileDirs = []
    # We build an array to pass into cat().
    for i in atoms:
        fileDirs.append(potcarDir + i + '/' + POTCAR_NAME)
    # print(fileDirs)

    cat(fileDirs, dirName, POTCAR_NAME)
    print('Successfully built POTCAR. Stored in {}'.format(dirName))

    return 0

# Build a POTCAR object using pymatgen, not testable unless on Odyssey
# TODO: test this on Odyssey.
def buildPotcar(poscarObj=getPoscarObj(), dirName=ROOT, useGivenPotcar=False): # Just specify useGivenPotcar and dirName (if not ROOT) if POTCAR file given.
    if useGivenPotcar:
        try:
            potcar = Potcar.from_file(checkPath(dirName) + POTCAR_NAME)
            return potcar
        except FileNotFoundError as err:
            print('Error:', err)
            print('Possible source of error: did you specify you want to use your own POTCAR without actually having one in the directory?')
            sys.exit()
        except Exception as err:
            print('Error:', err)
            print('Possible source of error: is your POTCAR file actually a valid POTCAR format?')
            sys.exit()

    else:
        try:
            atoms = poscarObj.site_symbols # Get array of atoms for potential, in order as given in POSCAR
        except Exception as err:
            print('Error:', err)
            print('Possible source of error: your POSCAR is likely wrong. You cannot import POTCAR until POSCAR is made properly for the PUC. Check POSCAR input and POTCAR maker function input.')
            sys.exit()
        #print(atoms)

        # Get a pymatgen Potcar object
        try:
            potcar = Potcar(atoms)
            print('POTCARs fetched by pymatgen. Writing to file...')
            potcar.write_file(checkPath(dirName) + POTCAR_NAME)
            return potcar
        except ValueError as err:
            print('Error:', err)
            print('\nPossible solution: run the following command as a one-time intialization to set up directories, if using Odyssey. \n\n\t{}\n'.format(POT_PMG_INIT_CMD))
            sys.exit()