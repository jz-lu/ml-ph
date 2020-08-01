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

