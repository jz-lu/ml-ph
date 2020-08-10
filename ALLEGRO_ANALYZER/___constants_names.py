# Program details

PRGM_COOL_NAME = '  __   __    __    ____  ___  ____   __         \n / _\\ (  )  (  )  (  __)/ __)(  _ \\ /  \\        \n/    \\/ (_/\\/ (_/\\ ) _)( (_ \\ )   /(  O )       \n\\_/\\_/\\____/\\____/(____)\\___/(__\\_) \\__/        \n  __   __ _   __   __    _  _  ____  ____  ____ \n / _\\ (  ( \\ / _\\ (  )  ( \\/ )(__  )(  __)(  _ \\\n/    \\/    //    \\/ (_/\\ )  /  / _/  ) _)  )   /\n\\_/\\_/\\_)__)\\_/\\_/\\____/(__/  (____)(____)(__\\_)'
PRGM_VERSION = 1.0
PRGM_END_CARD = ' ____  __ _  ____ \n(  __)(  ( \\(    \\\n ) _) /    / ) D (\n(____)\\_)__)(____/'

# Name of directory where this program is stored
# NOTE: this should not change! A command line argument specifies where the root directory of the calculation is, so we never need to move this
THIS_DIR = '/n/home04/jzlu/codes/ml-ph/ALLEGRO_ANALYZER/'

# Command line arguments in relaxation postprocessing
ELEDOS = 'eledos'
ELEBAND = 'eleband'
PHDOS = 'phdos'
PHBAND = 'phband'
ENERGIES = 'energies'
CMD_LINE_ARG_LIST = [ELEDOS, ELEBAND, PHDOS, PHBAND, ENERGIES]

# Directory names
START_BATCH_NAME = 'EXECUTABLE_BAT_DNE'
START_BATCH_PATH = '/n/home04/jzlu/codes/ml-ph/start_job_batch/' + START_BATCH_NAME
BATCH_FILE_NAME = 'STATIC_BAT_DNE'
BATCH_FILE_PATH = '/n/home04/jzlu/codes/ml-ph/ALLEGRO_ANALYZER/' + BATCH_FILE_NAME
RELAXATION_DIR_NAME = 'relaxation_calculations'
ANALYSIS_DIR_NAME = 'analyses'
TOTAL_ENER_DIR_NAME = 'energy_data'
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

KPOINTS_MESH_NAME = 'KPOINTS'
KPOINTS_LINE_NAME = 'LINE_KPOINTS'

VASP_RUN_XML_NAME = 'vasprun.xml'
VASP_RUN_OUT_NAME = 'vasp_%s.out'
VASP_RUN_ERR_NAME = 'vasp_%s.err'
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

