# Any miscellaneous constants and hard strings go here

unique_str = 'aUniqueStringThatIsVeryUnlikelyToBeReproducedRandomly_q345234572345'

# Error Messages
GENERAL_ERR_USAGE_MSG = 'Usage: python3 main.py <I/O DIRECTORY> <vdW T/F> <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband, energies).'
ERR_BAD_INPUT_DIR = 'Root I/O directory invalid. Check first command-line argument.'
ERR_BAD_DIR = 'Invalid directories.'
ERR_INVALID_VDW_FLAG = 'Invalid specification for van der Waals forces. Specify T or F in command line args.'
ERR_NO_POSCAR = 'Error: POSCAR required for input either invalid or not found.'
ERR_NO_POTCAR = 'No potentials found for atom given. Check potential list directory and verify POSCAR is correct.'
ERR_VASP_RUN_RELAX = 'Error running Vasp in relaxation:'
BAD_INPUT_ERR_MSG = 'Error: invalid command line arguments.'
ERR_BAD_KPOINTS_MODIFY_INPUT = 'Error: to modify KPOINTS file or Kpoints object you must choose samplingType as line or mesh.'
ERR_ENER_WRITE_FAIL = 'Failed to write out the energies to the specified root directory.'