# Name of directory where this program is stored
# NOTE: this should not change! A command line argument specifies where the root directory of the calculation is, so we never need to move this
THIS_DIR = '/Users/jonathanlu/Documents/ml-ph/ph-practice/gr-full'

# Command line arguments in relaxation postprocessing
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
ENERGIES = 'energies'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND, ENERGIES]

# Directory names
START_BATCH_NAME = 'EXECUTABLE_BAT_DNE'
START_BATCH_PATH = '/n/home04/jzlu/codes/ml-ph/static_batch_file/' + START_BATCH_NAME
BATCH_FILE_NAME = 'STATIC_BAT_DNE'
BATCH_FILE_PATH = '/n/home04/jzlu/codes/ml-ph/HTP/' + BATCH_FILE_NAME
PHONOPY_DIR_NAME = 'phonon_calculations'
PHDISP_STATIC_NAME = 'disp-'
PHDISP_DIR_NAME = PHDISP_STATIC_NAME + '%s/'

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
VASP_RUN_RELAX_OUT_NAME = 'vasp_relaxation.out'
VASP_RUN_RELAX_ERR_NAME = 'vasp_relaxation.err'
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

# Output names
ELEDOS_RAW_DATA_NAME = 'total_ele_dos_pmg'
ELEBAND_RAW_DATA_NAME = 'raw_data_ele_band_pmg'

# File formats
PLOT_FILE_FORMAT = 'eps'

