# All the constants and hard strings for vasp are here.
# NOTE: all directory strings must end in '/'

# Error checker for vasp.out
VASP_OUTFILE_LEN_THRESHOLD = 1000

### PARAMETERS FOR VASP INPUT FILES ###
## INCAR
INCAR_DEFAULT_SETTINGS = {'ISTART': 0, 
                          'ISMEAR': 0, 
                          'SIGMA': 0.01, 
                          'ENCUT': 800, 
                          'AMIN': 0.01, 
                          'NSW': 60, 
                          'EDIFF': 1e-08, 
                          'EDIFFG': -1e-08, 
                          'IBRION': 2, 
                          'SYMPREC': 0.0001, 
                          'NPAR': 6, 
                          'ADDGRID': True, 
                          'LREAL': False, 
                          'LWAVE': False, 
                          'LCHARG': False, 
                          'ALGO': 'N', 
                          'NELMIN': 5, 
                          'PREC': 'Accurate'} # We'll add in SYSTEM manually, EDIFFG can be larger if not phonon calc.
INCAR_VDW_SETTINGS = {'METAGGA': 'SCAN',
                      'LASPH': '.TRUE.',
                      'LUSE_VDW': '.TRUE.',
                      'BPARAM': 15.7}
#INCAR_DEFAULT_KEYS = ['SYSTEM', 'ISTART', 'ISMEAR', 'SIGMA', 'ENCUT', 'AMIN', 'NSW', 'EDIFF', 'EDIFFG', 'IBRION', 'SYMPREC', 'NPAR', 'ADDGRID', 'LREAL', 'LWAVE', 'LCHARG', 'ALGO', 'NELMIN', 'PREC']
#INCAR_VDW_KEYS = ['METAGGA', 'LASPH', 'LUSE_VDW', 'BPARAM']
ISTART = 0 # No WAVECAR input
ISMEAR = 0 # Gaussian smearing
SIGMA = 0.01 # Width of smear
ENCUT = 800 # Plane wave expansion cutoff index
AMIN = 0.01 # Parameter in initial approximation of dielectric function for screening
NSW = {'relax': 80, 'no_relax': 0, 'relax_low': 40, 'relax_very_low': 20} # Number of ionic relaxation steps
EDIFF = 1E-8 # Acceptable self-consistent energy difference in electronic relaxation 
EDIFFG = -1E-8 # Acceptable self-consistent energy difference (or forces if with negative sign as flag) in ionic relaxation
NEDOS = 6001 # Number of samples for electronic DOS calculation in the fixed energy range
IBRION = {'relax': 2, 'no_relax': -1} # 2 is a flag for do ionic relaxation, set to -1 when relaxation calculation is complete
ICHARG = {'default': 2, 'no_relax_sc': 1, 'no_relax_nsc': 11} # 2 is default initial charge density, 1 is self-consistent CONTCAR import, 11 is non self-consistent (fixed density) CONTCAR
SYMPREC = 0.0001
NPAR = 6
ADDGRID = '.TRUE.'
LREAL = '.FALSE.'
LWAVE = '.FALSE.'
LCHARG = '.FALSE.'
ALGO = 'N'
NELMIN = 5
PREC = 'Accurate'
# vdW settings
METAGGA = 'SCAN'
LASPH = '.TRUE.'
LUSE_VDW = '.TRUE.'
BPARAM = 15.7

## POTCAR
POT_PMG_INIT_CMD = 'cat “PMG_VASP_PSP_DIR: /n/kaxiras_lab/atomate_PPs” > ~/.pmgrc.yaml'
POT_NOPMG_DIR = '/n/kaxiras_lab/vasp.5.4.4/PPs/potpaw_PBE.54/' # This is the method of direct POTCAR without pymatgen!

## KPOINTS
RELAXATION_GRID_DENSITY = (21, 21, 1)
RELAXATION_GRID_SHIFT = (0, 0, 0)
NONRELAXATION_GRID_DENSITY = (21, 21, 1) # For precise calculations like DOS
NONRELAXATION_GRID_SHIFT = (0, 0, 0)
KPOINTS_LINE_INTS = 50 # Number of sampling k-points on each line
KPOINTS_LINE_HEX_STR = 'k-points along high symmetry lines\n50\nLine_mode\nReciprocal\n0.0 0.0 0.0 ! \Gamma\n0.5 0.5 0.0 ! M\n\n0.5 0.5 0.0 ! M\n0.666667 0.333333 0.0 ! K\n\n0.666667 0.333333 0.0 ! K\n0.0 0.0 0.0 ! \Gamma'

## POSCAR
FULL_RELAX_SELECTIVE_DYNAMICS_ARR = [True, True, True] # This is per-atom so you append this for every atom.
LAYER_RELAX_SELECTIVE_DYNAMICS_ARR = [False, False, True] # Same deal as above, only fixes interlayer spacing.