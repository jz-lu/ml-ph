# Command line arguments in relaxation postprocessing
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
ENERGIES = 'energies'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND, ENERGIES]

# Directories
## NOTE: it is crucial that all directories end in a forward slash. ##
ROOT = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full/' # File path to the scripts
BATCH_FILE_PATH = '/n/home04/jzlu/codes/ml-ph/STATIC_BAT_DNE'
DIR_ELEDOS = ROOT + ELEDOS + '/'
DIR_ELEBAND = ROOT + ELEBAND + '/'
PHONOPY_DIR_NAME = 'phonon_calculations'
PHDISP_STATIC_NAME = 'disp-'
PHDISP_DIR_NAME = PHDISP_STATIC_NAME + '%s/'
DIR_PHONOPY = ROOT + PHONOPY_DIR_NAME + '/'
DIR_PHDOS = DIR_PHONOPY + PHDOS + '/'
DIR_PHBAND = DIR_PHONOPY + PHBAND + '/'

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

VASP_RUN_XML_NAME = 'vasprun.xml'

