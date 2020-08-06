from __build_inputs import buildInitialInputs
from pymatgen.io.vasp.inputs import VaspInput

def getVaspObject(incar, kpoints, poscar, potcar, optionalFiles):
    return VaspInput(incar, kpoints, poscar, potcar, optional_files=optionalFiles)