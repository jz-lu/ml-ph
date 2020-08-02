from ___constants_names import POSCAR_UNIT_RELAXATION_NAME, ROOT
from ___constants_misc import ERR_NO_POSCAR
from __directory_searchers import find
from pymatgen.io.vasp.inputs import Poscar
import sys

# Get a POSCAR object
def getPoscarObj(poscarName=POSCAR_UNIT_RELAXATION_NAME, dirName=ROOT):
    poscarExists = find(poscarName, dirName)
    if poscarExists:
        return Poscar.from_file(dirName + poscarName)
    else:
        sys.exit(ERR_NO_POSCAR)
        return

# Get the total number of atoms in PUC
def getNumAtoms(poscarObj=getPoscarObj()):
    atomCountArr = poscarObj.natoms
    numAtoms = 0
    for i in atomCountArr:
        numAtoms += i
    return numAtoms

# Helps us label the input files
def getInputName(poscarObj=getPoscarObj()):
    sysName = ''
    # Labels it by atoms and number of each,  e.g. C-1, Si-2
    for i in range(0, len(poscarObj.site_symbols)):
        sysName = sysName + poscarObj.site_symbols[i] + '-' + str(poscarObj.natoms[i]) + ', '
    
    sysName = sysName[:-2] # trim off that last ', '

    return sysName