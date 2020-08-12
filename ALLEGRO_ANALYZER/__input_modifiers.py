from pymatgen.io.vasp.inputs import Incar, Kpoints, Poscar
import numpy as np

from ____exit_with_error import exit_with_error

from __directory_searchers import checkPath
from __query_inputs import getInputName, getNumAtoms

from ___constants_vasp import *
from ___constants_names import *
from ___constants_misc import ERR_BAD_KPOINTS_MODIFY_INPUT
from pymatgen.core.structure import Structure

# Functions for editing the input files depending on the calculation we need.
def modifyIncar(incar, addArr=None, delArr=None): # Add any elements and remove any elements respectively.
    # Function takes an incar object from pmg as input and we modify it.
    # Format of addArr: each element is a tuple (key, value) to add. Format of delArr: each element is just a key to delete.
    if addArr != None:
        for i in addArr:
            incar[i[0]] = i[1]
    
    if delArr != None:
        for i in delArr:
            if i in list(incar.keys()):
                incar.pop(delArr)
            else:
                print('Attempted deletion of %s failed as incar did not have that entry to begin with.'%(i))
    return incar

# Give the right Kpoints object for a different type of calculation, usually a denser sampling for DOS.
def modifyKpointsMesh(kpts_obj, meshDensity=NONRELAXATION_GRID_DENSITY, totShift=NONRELAXATION_GRID_SHIFT): # material name
    kpts_obj.kpts[0] = meshDensity
    kpts_obj.kpts_shift = totShift
    return kpts_obj

# This allows us to change the selective dynamics from relaxing only the interlayer spacing to full relaxation or vice versa.
def modifySelectiveDynamics(poscarObj, outDir, outFileName=POSCAR_UNIT_RELAXATION_NAME, relax_z_only=False, writeOut=False): # change writeOut and outDir when modifying in subfolders
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

# Perform a random lattice perturbation. If two dimensional, then we will need to do more work.
# TODO: check this function and how to modify it properly.
def randomAtomicPerturbations(distance, poscarObj, outDir, minDistance=None, twoDimensionalMat=True, writeOut=False, outFileName=POSCAR_UNIT_RELAXATION_NAME):
    # Random perturbation on sites in PUC. Choose a distance, and if fixed distance of perturbation do nothing else
    # If variable perturbation distance, set minDistance to something smaller than distance and we uniformly choose minDistance <= d <= distance
    
    # If the material only has 2-periodicity, then we can't allow the 3-oerturbation. 
    # So we store the projection of each unperturbed cooridnate onto 3rd lattice vector and change it back after pmg random perturbation.
    threeCoords = [] # Ordered Projections of atomic positions onto of 3rd basis vector, which is just scaled \hat{z}
    if twoDimensionalMat:
        for i in range(0, getNumAtoms(poscarObj)):
            threeCoords.append(poscarObj.structure.sites[i].frac_coords[2]) # last coordinate in frac_coords is proj on 3rd lattice vector

    poscarObj.structure.perturb(distance, min_distance=minDistance)
    print('Random perturbation of POSCAR-given sites complete. New object Poscar returned.')
    # If the numbers become huge, that is because the direction of perturbation pushed it into a neighbor, modulo back to the end of the PUC

    # Time to restore the 2D material's projection onto 3rd lattice vector to unperturbed state.
    if twoDimensionalMat:
        print('Perturbation is on a 2D solid. Restoring unperturbed state on 3rd dimension.')
        for i in range(0, len(threeCoords)):
            poscarObj.structure.sites[i].frac_coords[2] = threeCoords[i]

    if writeOut:
        print('Writing perturbed POSCAR to file...')
        poscarObj.write_file(checkPath(outDir) + outFileName)

# Make additions or deletions from the POSCAR PUC.
def addAtomToSite(atomName, dir_coords, poscarObj, outDir, writeOut=False, outFileName=POSCAR_UNIT_RELAXATION_NAME):
   try:
       newPoscarObj = poscarObj.structure.append(atomName, dir_coords, validate_proximity=True)
   except ValueError as err:
        print('Error:', err)
        exit_with_error('Suggested source of problem: you likely placed an atom that is too close to another one already in the PUC, which violates the obvious physics.')
   
   if writeOut:
       print('Writing updated POSCAR with inserted atom {} to file...'.format(atomName))
       newPoscarObj.write_file(checkPath(outDir) + outFileName)
   
   return newPoscarObj

# NOTE: this function requires the coordinates to be EXACT with the POSCAR i.e. same digit accuracy. e.g. 0.3333 will fail if POSCAR says 0.333333 or 0.33
def removeAtomFromSite(atomName, dir_coords, poscarObj, outDir, writeOut=False, outFileName=POSCAR_UNIT_RELAXATION_NAME):
    # Technically we could do it directly by index but in that case we may as well just use the direct pmg command
    # poscarObj.structure.remove_sites(array of indices here (start: 0) in the order of the POSCAR)

    # So instead we will let user input the atom and its position to be removed.
    index = None
    try:
        for i in range(0, getNumAtoms(poscarObj)):
            if np.array_equal(poscarObj.structure.sites[i].frac_coords, dir_coords): # no issue
                index = i
        if index == None:
            exit_with_error('Error: the atom %s at (%f, %f, %f) specified for removal does not exist in the Poscar object.'%(atomName, dir_coords[0], dir_coords[1], dir_coords[2]))

        newPoscarObj = poscarObj.structure.remove_sites([index])
        print('Removed atom %s at coordinates (%f, %f, %f) from Poscar object.'%(atomName, dir_coords[0], dir_coords[1], dir_coords[2]))
    except Exception as err:
        print('Error:', err)
        
    print(poscarObj.structure.sites)

    if writeOut:
       print('Writing updated POSCAR with removed atom %s at (%f, %f, %f) to file...'%(atomName, dir_coords[0], dir_coords[1], dir_coords[2]))
       newPoscarObj.write_file(checkPath(outDir) + outFileName)
    
    return newPoscarObj

# The following two functions prepare incar for nonrelaxed calculations, respectively self-consistent and nonself-consistent for dos/ph and band.
def getSelfConNoRelIncar(incarObj):
    settings_to_add = [('ICHARG', ICHARG['no_relax_sc']), 
                       ('IBRION', IBRION['no_relax']), 
                       ('NSW', NSW['no_relax']),
                       ('SIGMA', SIGMA['wide'])]
    incarObj = modifyIncar(incarObj, settings_to_add)
    return incarObj

def getNonSelfConNoRelIncar(incarObj):
    settings_to_add = [('ICHARG', ICHARG['no_relax_nsc']),
                       ('IBRION', IBRION['no_relax']), 
                       ('NSW', NSW['no_relax']),
                       ('LCHARG', LCHARG['no_write_charge'])] # No NEDOS, thats just for plotting DOS
    incarObj = modifyIncar(incarObj, settings_to_add)
    return incarObj