# Any miscellaneous constants and hard strings go here

unique_str = 'aUniqueStringThatIsVeryUnlikelyToBeReproducedRandomly_q345234572345'

# Error Messages
ERR_BAD_DIR = 'Invalid directories.'
ERR_NO_POSCAR = 'Error: no input POSCAR found.'
ERR_NO_POTCAR = 'No potentials found for atom given. Check potential list directory and verify POSCAR is correct.'
GENERAL_ERR_USAGE_MSG = 'Usage: python3 main.py <I/O DIRECTORY> <arg1> <arg2> .... Specify at least one arg (eledos, eleband, phdos, phband, energies).'
BAD_INPUT_ERR_MSG = 'Error: invalid command line arguments.'
ERR_BAD_KPOINTS_MODIFY_INPUT = 'Error: to modify KPOINTS file or Kpoints object you must choose samplingType as line or mesh.'