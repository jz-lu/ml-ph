# Command line arguments in relaxation postprocessing
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND]

# Directories
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full/' # File path to the scripts
## Note!! It is crucial that ROOT end in a forward slash. ##
DIR_ELEDOS = ROOT + ELEDOS + '/'
DIR_ELEBAND = ROOT + ELEBAND + '/'
DIR_PHDOS = ROOT + PHDOS + '/'
DIR_PHBAND = ROOT + PHBAND + '/'

# File names
POSCAR_NAME = 'POSCAR'
KPOINTS_NAME = 'KPOINTS'
INCAR_RELAXATION_NAME = 'INCAR'
POTCAR_NAME = 'POTCAR'
CONTCAR_NAME = 'CONTCAR'
CHGCAR_NAME = 'CHGCAR'

POSCAR_UNIT_RELAXATION_NAME = 'POSCAR'
POSCAR_UNIT_PHONOPY_NAME = 'POSCAR_unit'

KPOINTS_MESH_NAME = 'KPOINTS'
KPOINTS_LINE_NAME = 'LINE_KPOINTS'


