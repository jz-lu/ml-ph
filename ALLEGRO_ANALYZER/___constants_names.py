# Program details

PRGM_COOL_NAME = '\n█▀▄▀█ █░░ ▄▄ █▀█ █░█ ▀   █▀ ▀█▀ ▄▀█ █▀█ ▀█▀ █ █▄░█ █▀▀\n█░▀░█ █▄▄ ░░ █▀▀ █▀█ ▄   ▄█ ░█░ █▀█ █▀▄ ░█░ █ █░▀█ █▄█\n'
PRGM_VERSION = 1.0
PRGM_END_CARD = '\n\n█▀▄▀█ █░░ ▄▄ █▀█ █░█   █▀▀ ▀▄▀ █ ▀█▀ █ █▄░█ █▀▀\n█░▀░█ █▄▄ ░░ █▀▀ █▀█   ██▄ █░█ █ ░█░ █ █░▀█ █▄█\n\n'
BACKUP_END_CARD = ' ____  __ _  ____ \n(  __)(  ( \\(    \\\n ) _) /    / ) D (\n(____)\\_)__)(____/'

# Name of directory where this program is stored
# NOTE: this should not change! A command line argument specifies where the root directory of the calculation is, so we never need to move this
THIS_DIR = '/n/home04/jzlu/codes/ml-ph/start_job_batch/'

# Default input specifications filename
DEFAULT_INPUT_FILENAME = '/n/home04/jzlu/codes/ml-ph/inputfiles/input.txt'
DIR_SHORTCUTS = ['.', './', '..', '../']

# Command line arguments in relaxation postprocessing
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
ENERGIES = 'energies'
PH = 'ph'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND, ENERGIES, PH]

# Calculation types
TYPE_RELAX_BASIC = 0
TYPE_RELAX_CONFIG = 1
TYPE_TWISTED_CONFIG = 2
TYPE_NORELAX_BASIC = 3
TYPE_FLAGS = [TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC]

# Directory and batch names
START_BATCH_NAME = 'EXECUTABLE_BAT_DNE'
CONFIG_BATCH_NAME = 'SHIFT_BAT_DNE_'
START_BATCH_PATH = '/n/home04/jzlu/codes/ml-ph/start_job_batch/' + START_BATCH_NAME
START_BATCH_OUTFILE_START = 'job_'
BATCH_FILE_NAME = 'STATIC_BAT_DNE'
CODE_DIR = '/n/home04/jzlu/codes/ml-ph/ALLEGRO_ANALYZER/'
BATCH_FILE_PATH = CODE_DIR + BATCH_FILE_NAME
RELAXATION_DIR_NAME = 'relaxation'
ANALYSIS_DIR_NAME = 'analyses'
TOTAL_ENER_DIR_NAME = 'energy_data'
PHONOPY_DIR_NAME = 'phonon'
PHDISP_STATIC_NAME = 'disp'
PHDISP_DIR_NAME = PHDISP_STATIC_NAME + '%d/'
MY_EMAIL = 'jlu@college.harvard.edu'
CONFIG_DATA_DIR = 'config_data/'
CONFIG_SUBDIR_NAME = 'shift_'
SHIFT_NAME = 'shiftcoord.txt'
COB_NPY_NAME = 'cob.npy' # must end in .npy
LIDXS_NPY_NAME = 'lidxs.npy' # must end in .npy
DSAMPLE_ENERGIES_PRE = 'd_energies'
DSAMPLE_ENERGIES_NAME = DSAMPLE_ENERGIES_PRE + '.npy'
DSAMPLE_ENERGIES_TXT = DSAMPLE_ENERGIES_PRE + '.txt'
DSAMPLE_SPACINGS_PRE = 'd_spacings'
DSAMPLE_SPACINGS_NAME = DSAMPLE_SPACINGS_PRE + '.npy'
DSAMPLE_SPACINGS_TXT = DSAMPLE_SPACINGS_PRE + '.txt'
DSAMPLE_FORCES_PRE = 'd_forces'
DSAMPLE_FORCES_NAME = DSAMPLE_FORCES_PRE + '.npz'
FGSFE_COEFF_PRE = 'gsfe_coef'
FGSFE_COEFF_NAME = FGSFE_COEFF_PRE + '.npy'
FGSFE_COEFF_TXT = FGSFE_COEFF_PRE + '.txt'
FGSFE_SCORE_NAME = 'gsfe_score.npy'
FGSFE_PLOT_NAME = 'gsfe_pva.png'
MONOLAYER_DIR_NAME = 'layer_'
CONFIG_DIR_NAME = 'config'

# File names
POSCAR_NAME = 'POSCAR'
SPOSCAR_NAME = 'SPOSCAR'
KPOINTS_NAME = 'KPOINTS'
INCAR_RELAXATION_NAME = 'INCAR'
POTCAR_NAME = 'POTCAR'
CONTCAR_NAME = 'CONTCAR'
OSZICAR_NAME = 'OSZICAR'
CHGCAR_NAME = 'CHGCAR'

POSCAR_UNIT_RELAXATION_NAME = 'POSCAR'
POSCAR_CONFIG_NAMEPRE = 'POSCAR_'
POSCAR_PH_NAMEPRE = 'POSCAR-'

KPOINTS_MESH_NAME = 'KPOINTS'
KPOINTS_LINE_NAME = 'LINE_KPOINTS'

VASP_RUN_XML_NAME = 'vasprun.xml'
VASP_RUN_OUT_NAME = '%s.out'
VASP_RUN_ERR_NAME = '%s.err'
DOSCAR_NAME = 'DOSCAR'
OUTCAR_NAME = 'OUTCAR'

OUTPUT_DIR_NAME = 'outputs'
COMBINED_ELE_OUTPUTS_NAME = 'additional_outputs'

TOT_ENERGIES_NAME = 'total_energies.txt'

PH_FORCE_SETS_NAME = 'FORCE_SETS'
PH_BASIC_INFO_FILENAME = 'phonopy_basic_structure_input.txt'
PH_MESH_CONF_NAME = 'mesh.conf'
PH_BAND_CONF_NAME = 'band.conf'
PH_BAND_RAW_DATA_NAME = 'ph_dos_raw'
PH_DISP_YAML_NAME = 'phonopy_disp.yaml'
DEFAULT_PH_BAND_PLOT_NAME = 'phband.png'

# Output names
ELEDOS_RAW_DATA_NAME = 'total_ele_dos_pmg'
ELEBAND_RAW_DATA_NAME = 'raw_data_ele_band_pmg'

# File formats
PLOT_FILE_FORMAT = 'eps'

