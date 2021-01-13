# Any miscellaneous constants and hard strings go here

# Number of available cores for muiltiprocessing
NUM_AVAILABLE_CORES = 12

unique_str = 'aUniqueStringThatIsVeryUnlikelyToBeReproducedRandomly_q345234572345'

# Any configuration space errors
ERR_INCONSISTENT_NORMS = 'Inconsistent norms of in-plane lattice vectors for one of the input POSCARS'
ERR_INCONSISTENT_LATTICES = 'Error: set of lattice vectors are not identical to the first lattice vector, up to lattice constant.'
ERR_ATOMS_OUT_OF_PLANE = 'Error: all input POSCARs must have atoms zeroed wrt the out-of-plane (3rd basis) vector. \n\n\tFound problem (first occurrence) in POSCAR file %d on atom %d.\n\n'
ERR_WRONG_ZVEC = 'Error: at least one lattice z-vector (i.e. 3rd basis vector) is inconsistent with constant fixed length %f.'
ERR_INVALID_GRID = 'Error: invalid grid. Usage: (a, b, 1). Note 1 since only 2D materials are supported.'

# Error Messages
GENERAL_ERR_USAGE_MSG = 'Usage: python3 start.py <typeFlag: 0 for single, 1 for config> <I/O DIRECTORY> <vdW T/F> <GAMMA or MP> <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband, energies).\n\n\t If using input file...\t Usage: "python3 start.py -f <filename>" where parameters are in the input file separated by a newline\n\n'
ERR_BAD_TYPE_FLAG = 'Type flag not correctly specified. Select 0 for normal calculation and 1 for configuration space sampling.'
ERR_BAD_INTERLAYER_DISTANCE = 'Interlayer distance invalid.'
ERR_BAD_INPUT_DIR = 'Unable to find given root I/O directory: directory must exist and have input CAR files.'
ERR_BAD_DIR = 'Invalid directories.'
ERR_INVALID_VDW_FLAG = 'Invalid specification for van der Waals forces. Specify T or F in command line args.'
ERR_NO_INPUT = 'Input not imported properly to handler class. No initialization done in instance of InputData class. Check code workflow.'
ERR_NO_POSCAR = 'Error: POSCAR required for input either invalid or not found.'
ERR_NO_POTCAR = 'No potentials found for atom given. Check potential list directory and verify POSCAR is correct.'
ERR_INVALID_FINDDIR_PARAM = 'Error: invalid specification for type of search (searchType) in calling filesInDir.'

ERR_VASP_RUN_RELAX = 'Error running Vasp in relaxation:'
ERR_VASP_NOT_CONVERGED = 'Error in running Vasp relaxation. Too many attempts made to converge relaxation, all failed. Check results of calculation for details.'
BAD_INPUT_ERR_MSG = 'Error: invalid command line arguments.'
ERR_BAD_KPOINTS_MODIFY_INPUT = 'Error: to modify KPOINTS file or Kpoints object you must choose samplingType as line or mesh.'
ERR_ENER_WRITE_FAIL = 'Failed to write out the energies to the specified root directory.'

ERR_PH_FORCE_SETS = 'An error occurred while generating FORCE_SETS with phonopy.'
ERR_PH_FORCE_SETS_NOT_FOUND = 'FORCE_SETS not found in the directory. Check logs for the creation of FORCE_SETS by phonopy. The function in this script is a good place to start.'

ERR_PH_CANNOT_GEN_MESHCONF = 'Error in creating mesh.conf. Invalid KPOINTS object used. Did you use a line KPOINTS instead of mesh?'