from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
from __dirModifications import checkPath
from __query_inputs import *

from ___constants_vasp import FULL_RELAX_SELECTIVE_DYNAMICS_ARR, LAYER_RELAX_SELECTIVE_DYNAMICS_ARR, KPOINTS_MESH_DOS
from ___constants_names import *

# Functions for editing the input files depending on the calculation we need.
def modifyIncar(incar, addArr, delArr): # Add any elements and remove any elements respectively.
    # Function takes an incar object from pmg as input and we modify it.
    # Format of addArr: each element is a tuple (key, value) to add. Format of delArry: each element is just a key to delete.
    for i in addArr:
        incar[i[0]] = i[1]
    for i in delArr:
        incar.pop(delArr)
    return incar

# Give the right Kpoints object so that we can write it into the proper subdirectory.
def modifyKpoints(samplingType, matName, totShift=(0, 0, 0)): # material name
    if samplingType == 'mesh':
        kpoints_new = Kpoints.gamma_automatic( KPOINTS_MESH_DOS, totShift )
        kpoints_new.comment = 'Kpoints grid for ' + getInputName()
    else:
        # Generate a line KPOINTS for band calculations. See documentation of the library for details.
        # TODO: for now, we just import a LINE_KPOINTS file since the autogenerator doesn't seem to work.
        kpoints_new = Kpoints.from_file(ROOT + '/LINE_KPOINTS')
    return kpoints_new

# This allows us to change the selective dynamics from relaxing only the interlayer spacing to full relaxation or vice versa.
def modifySelectiveDynamics(poscarObj=getPoscarObj(), relax_z_only=False, writeOut=False, outDir=ROOT, outFileName=POSCAR_UNIT_RELAXATION_NAME): # change writeOut and outDir when modifying in subfolders
    natoms = getNumAtoms(poscarObj)
    sdMatrix = [] # Nx3 array of booleans for selective dynamic flags

    if relax_z_only:
        for i in range(0, natoms):
            sdMatrix.append(LAYER_RELAX_SELECTIVE_DYNAMICS_ARR)
        poscarObj.selective_dynamics = sdMatrix
        print('Updated selective dynamics of POSCAR.')
    else:
        if type(poscarObj.selective_dynamics) is list: # if it has no selective dynamics, then it will be NoneType not list and we don't need to do anything
            for i in range(0, natoms):
                sdMatrix.append(FULL_RELAX_SELECTIVE_DYNAMICS_ARR)
            poscarObj.selective_dynamics = sdMatrix
            print('Updated selective dynamics of POSCAR.')
        else:
            print('Selective dynamics for POSCAR already match desired settings. No updates performed.')
    
    if writeOut:
        # We can write it out into the outDir to update POSCAR
        print('Writing out POSCAR with desired selective dynamics settings...')
        poscarObj.write_file(checkPath(outDir) + outFileName)
    
    return poscarObj

def perturbAtoms():
    print('foo')

def addAtoms():
    print('bar')