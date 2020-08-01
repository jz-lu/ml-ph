from pymatgen.io.vasp.inputs import Incar, Kpoints

# Any relevant constants/hard strings
KPOINTS_MESH_DOS = (21, 21, 1) # Mesh sampling
KPOINTS_LINE_INTS = 50 # Number of sampling k-points on each line
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full' # File path to the scripts

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
        # TODO: add kpoints_new.comments = (species info + mesh)
    else:
        # Generate a line KPOINTS for band calculations. See documentation of the library for details.
        # TOCHECK: for now, we just import a LINE_KPOINTS file since the autogenerator doesn't seem to work.
        kpoints_new = Kpoints.from_file(ROOT + '/LINE_KPOINTS')
    return kpoints_new
    