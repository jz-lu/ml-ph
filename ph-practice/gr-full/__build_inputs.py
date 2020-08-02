from ___constants_misc import *
from ___constants_names import *
from ___constants_vasp import *
from __directory_searchers import find
from __dirModifications import move, cat
import sys

from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints # No Potcar class, we'll just do it with cat the ol' fashioned way

# Get a POSCAR object
def getPoscarObj(poscarName=POSCAR_UNIT_RELAXATION_NAME, dirName=ROOT):
    poscarExists = find(poscarName, dirName)
    if poscarExists:
        return Poscar.from_file(dirName + poscarName)
    else:
        sys.exit(ERR_NO_POSCAR)
        return

# Helps us label the input files
def getInputName(dirName=ROOT, poscarName=POSCAR_UNIT_RELAXATION_NAME):
    poscar = Poscar.from_file(ROOT + poscarName)
    sysName = ''
    for i in poscar.site_symbols:
        sysName = sysName + i + '-'
    sysName += str(poscar.natoms[0])

    return sysName

# Returns pymatgen Kpoints object
def buildRelaxationKpoints(dirName=ROOT, grid=RELAXATION_GRID_DENSITY, shift=RELAXATION_GRID_SHIFT, isRelaxed=False, writeOut=True):
    kpointsExists = find(KPOINTS_MESH_NAME, ROOT) # Look for existing KPOINTS
    if kpointsExists:
        kpoints = Kpoints.from_file(ROOT + KPOINTS_MESH_NAME)
        print('KPOINTS file inputted by user. Importing input for computation...')
    else:
        if grid == RELAXATION_GRID_DENSITY and shift == RELAXATION_GRID_SHIFT:
            print('No KPOINTS file found. Generating default relaxation mesh...')
        else:
            print('No KPOINTS file found. Generating mesh according to user inputted grid...')
        kpoints = Kpoints.gamma_automatic(grid, shift) # Only supports Gamma centered Monkhorst Pack for now
    
    kpoints.comment = getInputName() # Standardize labeling style

    if (not kpointsExists) and writeOut:
        print('Writing new KPOINTS file...')
        kpoints.write_file(ROOT + KPOINTS_MESH_NAME)
    
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
        incar = Incar.from_file(ROOT + INCAR_RELAXATION_NAME)
        for i in range(0, len(default_keys)):
            if default_keys[i] not in incar:
                print('{} not in INCAR'.format(i))
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

# Writes a properly constructed POTCAR (not an object)
# Returns 0 on success, 1 on failure
def buildPotcar_noPMG(dirName=ROOT, potcarDir=POT_NOPMG_DIR, poscarObj=getPoscarObj()):
    atoms = poscarObj.site_symbols
    numAtoms = len(atoms)

    if numAtoms == 1:
        # Just get the one POTCAR
        potcarExists = find(atoms[0], POT_NOPMG_DIR)
        if not potcarExists:
            sys.exit(ERR_NO_POTCAR)
        movePotcar = move('POTCAR', POT_NOPMG_DIR + atoms[0], ROOT)
        if movePotcar.stderr !=  '':
            print('Error in moving POTCARs to {}. Check process by hand.'.format(ROOT))
            return 1
    
    elif numAtoms > 1:
        # We need to concatenate POTCARs.First we do the first two then the rest.
        
        
    