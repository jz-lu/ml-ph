from ___constants_misc import *
from ___constants_names import *
from ___constants_vasp import *
from __directory_searchers import find

from pymatgen.io.vasp.inputs import Poscar, Incar, Kpoints # No Potcar class, we'll just do it with cat the ol' fashioned way

def getInputName(dirName=ROOT, poscarName=POSCAR_UNIT_RELAXATION_NAME):
    poscar = Poscar.from_file(ROOT + poscarName)
    sysName = ''
    for i in poscar.site_symbols:
        sysName = sysName + i + '-'
    print(poscar.natoms)
    sysName += str(poscar.natoms[0])

    return sysName

# Returns pymatgen Kpoints object
def buildRelaxationKpoints(dirName=ROOT, grid=RELAXATION_GRID_DENSITY, totShift=RELAXATION_GRID_SHIFT, isRelaxed=False, writeOut=True):
    kpointsExists = find(KPOINTS_MESH_NAME, ROOT) # Look for existing KPOINTS
    if kpointsExists:
        kpoints = Kpoints.from_file(ROOT + KPOINTS_MESH_NAME)
        print('KPOINTS file inputted by user. Importing input for computation...')
    else:
        if grid == RELAXATION_GRID_DENSITY and totShift == RELAXATION_GRID_SHIFT:
            print('No KPOINTS file found. Generating default relaxation mesh...')
        else:
            print('No KPOINTS file found. Generating mesh according to user inputted grid...')
        kpoints = Kpoints.gamma_automatic(grid, totShift) # Only supports Gamma centered Monkhorst Pack for now
    
    kpoints.comment = getInputName() # Standardize labeling style

    if (not kpointsExists) and writeOut:
        print('Writing new KPOINTS file...')
        kpoints.write_file(ROOT + KPOINTS_MESH_NAME)
    
    if not writeOut:
        print('KPOINTS object created without file writing.')

    return kpoints

buildRelaxationKpoints()

    
# Returns pymatgen Incar object
def buildRelaxationIncar(dirName=ROOT, vdW=True, writeOut=True):
    # First let's build an array of default keys, which is helpful later
    default_keys = list(INCAR_DEFAULT_SETTINGS.keys())
    vdW_keys = list(INCAR_VDW_SETTINGS.keys())

    incarExists = find(INCAR_RELAXATION_NAME, ROOT) # Look for existing INCAR
    if not incarExists:
        # Create a default
        incar = Incar.from_string('')

# incar = Incar.from_file('INCAR')
# print(incar.items())
# print(incar.values())

# printout = {}
# for i in range(0, len(INCAR_DEFAULT_KEYS)):
#     printout[INCAR_DEFAULT_KEYS[i]] = (['C:graphene-relax-phonon', 0, 0, 0.01, 800, 0.01, 1, 1e-08, -1e-08, -1, 0.0001, 6, True, False, False, False, 'N', 5, 'Accurate'])[i]
# print(printout)
# for i in range(0, len(INCAR_DEFAULT_SETTINGS)):
#     print(list(INCAR_DEFAULT_SETTINGS.keys())[i], '=', list(INCAR_DEFAULT_SETTINGS.values())[i])