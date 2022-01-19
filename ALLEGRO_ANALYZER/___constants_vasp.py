# All the constants and hard strings for vasp are here.
# NOTE: all directory strings must end in '/'
from ___constants_phonopy import SUPER_DIM
from math import ceil

# Unit conversion
VASP_FREQ_TO_INVCM_UNITS = 15.633302*33.356
USE_SCAN_VDW = True # change to false if using DFT-D3

# POSCAR z-basis vector magnitude, increase if we add a lot more (>3) layers
# Z_LATTICE_SIZE = 20.0
# Z_LAYER_SEP = 0.376723775 # Layer separation in terms of the z-basis vector
Z_LATTICE_SIZE = 35.0
Z_LAYER_SEP = {
  'WSe2': 0.2045,
  'MoSe2': 0.2045,
  'MoS2_0': 0.19405515, 
  'MoS2_180': 0.19403429,
  'Gr': 0.103
} # Layer separation in terms of the z-basis vector
Z_LAYER_SEP = Z_LAYER_SEP['MoSe2']
# Precision of POSCAR lattice and atoms, in number of decimal places
POSCAR_PRECISION = 6
POSCAR_PREC_COMP_THRESHOLD = 0.0001
INIT_EDIFF0 = 1e-3

# Error checker for vasp.out
VASP_OUTFILE_LEN_THRESHOLD = 1000
VASP_MAX_CONVERGENCE_ATTEMPTS = {
  'L': 6, 'T': 2
} # Number of times we rerun the relaxation before we quit and say we failed to relax.
VASP_OUTFILE_CONVERGENCE_LINE_CUT = 6 # Number of lines from the bottom we scan looking for an error message

### PARAMETERS FOR VASP INPUT FILES ###
## INCAR
ISTART = 0 # No WAVECAR input
ISMEAR = 0 # Gaussian smearing
SIGMA = {'wide': 0.1, 'narrow': 0.05} # Width of smear
ENCUT = 800 # Plane wave expansion cutoff index
AMIN = 0.01 # Parameter in initial approximation of dielectric function for screening
NSW = {'relax': 300, 'no_relax': 0, 'relax_low': 100, 'relax_very_low': 80} # Number of ionic relaxation steps
EDIFF = {'relax': 1E-8, 'vdw_relax': 1E-6, 'no_relax': 1E-8} # Acceptable self-consistent energy difference in electronic relaxation 
EDIFFG = {'relax': -1E-8, 'vdw_relax': -1E-6, 'no_relax': -1E-8} # Acceptable self-consistent RELAXATION energy difference (or forces if with negative sign as flag) in ionic relaxation
NEDOS = 6001 # Number of samples for electronic DOS calculation in the fixed energy range
IBRION = {'relax': 2, 'no_relax': -1, 'dfpt': 8} # set to -1 to indicate not to relax
ICHARG = {'default': 2, 'no_relax_sc': 1, 'no_relax_nsc': 11} # 2 is default initial charge density, 1 is self-consistent CONTCAR import, 11 is non self-consistent (fixed density) CONTCAR
SYMPREC = 0.0001
NPAR = 2
ADDGRID = '.TRUE.'
LREAL = '.FALSE.'
LWAVE = '.FALSE.'
LCHARG = {'write_charge': '.TRUE.', "no_write_charge": '.FALSE.'}
ALGO = 'N'
NELMIN = 5
PREC = 'Accurate'
POTIM = 0.1

# SCAN-rvv10 vdW settings
METAGGA = 'R2SCAN'
LASPH = '.TRUE.'
LUSE_VDW = '.TRUE.'
BPARAM = 15.7

# DFT-D3 vdW settings
GGA = 'PE'
IVDW = 12

INCAR_RELAX_SETTINGS = {'ISTART': ISTART, 
                          'ISMEAR': ISMEAR, 
                          'SIGMA': SIGMA['narrow'], 
                          'ENCUT': ENCUT, 
                          'AMIN': AMIN, 
                          'NSW': NSW['relax'], 
                          'EDIFF': EDIFF['relax'], 
                          'EDIFFG': EDIFFG['relax'], 
                          'IBRION': IBRION['relax'], 
                          'NPAR': NPAR, 
                          'ADDGRID': ADDGRID, 
                          'LREAL': LREAL, 
                          'LWAVE': LWAVE, 
                          'LCHARG': LCHARG['write_charge'], 
                          'ALGO': ALGO, 
                          'NELMIN': NELMIN, 
                          'PREC': PREC, 
                          'POTIM': POTIM
                        } # We'll add in SYSTEM manually, EDIFFG can be larger if not phonon calc.
INCAR_VDW_SETTINGS = {
    'GGA': GGA,
    'IVDW': IVDW
}
if USE_SCAN_VDW:
  INCAR_VDW_SETTINGS = {'METAGGA': METAGGA,
                        'LASPH': LASPH,
                        'LUSE_VDW': LUSE_VDW,
                        'BPARAM': BPARAM, 
                        'POTIM': 6*POTIM, 
                        'EDIFF': EDIFF['vdw_relax'], 
                        'EDIFFG': EDIFFG['vdw_relax']
                        }
INCAR_NORELAX_SCON_SETTINGS = {'ISTART': ISTART, 
                          'ISMEAR': ISMEAR, 
                          'SIGMA': SIGMA['narrow'], 
                          'ENCUT': ENCUT, 
                          'AMIN': AMIN, 
                          'NSW': NSW['no_relax'], 
                          'EDIFF': EDIFF['no_relax'], 
                          'EDIFFG': EDIFFG['no_relax'], 
                          'IBRION': IBRION['no_relax'], 
                          'NPAR': NPAR, 
                          'ADDGRID': ADDGRID, 
                          'LREAL': LREAL, 
                          'LWAVE': LWAVE, 
                          'LCHARG': LCHARG['no_write_charge'], 
                          'ALGO': ALGO, 
                          'NELMIN': NELMIN, 
                          'PREC': PREC}


## POTCAR
POT_PMG_INIT_CMD = 'cat “PMG_VASP_PSP_DIR: /n/kaxiras_lab/atomate_PPs” >> ~/.pmgrc.yaml'
POT_NOPMG_DIR = '/n/kaxiras_lab/vasp.5.4.4/PPs/potpaw_PBE.54/' # This is the method of direct POTCAR without pymatgen!

## KPOINTS
RELAXATION_GRID_DENSITY = (17, 17, 1)
RELAXATION_GRID_SHIFT = (0, 0, 0)
NONRELAXATION_GRID_DENSITY = (21, 21, 1)
SUPERCELL_GRID_DENSITY = (5, 5, 1)
NONRELAXATION_GRID_SHIFT = (0, 0, 0)
PHONOPY_GRID_DENSITY = (15, 15, 1)
PHONOPY_GRID_SHIFT = (0, 0, 0)
KPOINTS_LINE_INTS = 50 # Number of sampling k-points on each line
# pylint: disable=anomalous-backslash-in-string
KPOINTS_LINE_HEX_STR = 'k-points along high symmetry lines\n50\nLine_mode\nReciprocal\n0.0 0.0 0.0 ! \Gamma\n0.5 0.5 0.0 ! M\n\n0.5 0.5 0.0 ! M\n0.666667 0.333333 0.0 ! K\n\n0.666667 0.333333 0.0 ! K\n0.0 0.0 0.0 ! \Gamma'

## POSCAR
FULL_RELAX_SELECTIVE_DYNAMICS_ARR = [True, True, True] # Free movement
NO_RELAX_SELECTIVE_DYNAMICS_ARR = [False, False, False] # Fixes position
LAYER_RELAX_SELECTIVE_DYNAMICS_ARR = [False, False, True] # Fixes in-plane position
