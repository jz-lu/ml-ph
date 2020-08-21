# Functions needed: importing the poscar object, finding the first layer, sampling the grid, and 
# Exercise: write a class constructed on a POSCAR object (defaulted to taking from file) and have a bunch of member functions on it that do above

from ___constants_names import POSCAR_NAME
from __directory_searchers import checkPath
from pymatgen.io.vasp.inputs import Poscar

class ConfigSpaceSampling:
    def __init__(self, dirName, poscar=None):
        dirName = checkPath(dirName)
        if poscar == None:
            self.poscar = Poscar.from_file(dirName + POSCAR_NAME)
        else:
            self.poscar = poscar
        self.dirName = dirName

    