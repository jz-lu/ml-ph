# Constants related to phonopy

# Preprocessing
POSCAR_UNIT_NAME = 'POSCAR_unit' # Unit cell POSCAR name
PHONOPY_NUM_LINE_INTS = 101 # Number of samples along each line in kpoints path, should be an odd number so the center lines up
SUPER_DIM = (3, 3, 1)
LARGE_SUPER_DIM = (6,6,1)
SUPER_DIM_STR = ' '.join(map(str, SUPER_DIM)) # Supercell size
LARGE_SUPER_DIM_STR = ' '.join(map(str, SUPER_DIM)) # Supercell size
PHONOPY_DISP_MSG = 'Phonopy activated. Getting displacements...'
PHONOPY_ORG_MSG = 'Organizing displacement files into respective subdirectories...'
PHONOPY_DISP_ERR_1 = 'Error: failed to find any displacement files. Phonopy likely could not interpret input. Check input and try again.'
PHONOPY_DISP_ERR_2 = 'Error: unreasonably large number of displacement files (>999). Check input.'
PH_ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full' # path to the phonopy files
FC_CONF_NAME = 'fc.conf'