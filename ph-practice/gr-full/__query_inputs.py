from ___constants_misc import ERR_NO_POSCAR

from __directory_searchers import find, checkPath
from pymatgen.io.vasp.inputs import Poscar, Kpoints
import sys
    
# Get the total number of atoms in PUC
def getNumAtoms(poscarObj):
    atomCountArr = poscarObj.natoms
    numAtoms = 0
    for i in atomCountArr:
        numAtoms += i
    return numAtoms

# Helps us label the input files
def getInputName(poscarObj):
    sysName = ''
    # Labels it by atoms and number of each,  e.g. C-1, Si-2
    for i in range(0, len(poscarObj.site_symbols)):
        sysName = sysName + poscarObj.site_symbols[i] + '-' + str(poscarObj.natoms[i]) + ', '
    
    sysName = sysName[:-2] # trim off that last ', '

    return sysName

# Same thing but no spaces
def getInputFormula(poscarObj):
    sysName = ''
    # Labels it by atoms and number of each,  e.g. C-1, Si-2
    for i in range(0, len(poscarObj.site_symbols)):
        sysName = sysName + poscarObj.site_symbols[i] + str(poscarObj.natoms[i])

    return sysName