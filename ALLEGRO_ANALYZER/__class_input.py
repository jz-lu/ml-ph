from ____exit_with_error import exit_with_error
from ___constants_vasp import BPARAM
from ___constants_misc import *
from ___constants_config import GRID_SAMPLE_HIGH, GRID_SAMPLE_LOW, DIAG_SAMPLE_HIGH, DIAG_SAMPLE_LOW
from ___constants_names import (
    CMD_LINE_ARG_LIST, 
    ENERGIES, ELEBAND, PHBAND, PHDOS, PH, 
    TYPE_FLAGS, CFG_STRS, CFG_DIAG, CFG_Z, CFG_LC, 
    TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC, 
    ISPC_SAMPLE_INAME, LC_SAMPLE_INAME
)
from ___constants_vasp import EDIFF
from __directory_searchers import checkPath
import os, argparse
import numpy as np

# Since this is a simple parsing class, we'll dispense with the usually reasonable move of mangling member variables to simulate privacy.
class InputData:
    input_imported = False

    def __init__(self, args):
        self.super_dim = args.super; assert args.super > 0
        print(f"Supercell size: {self.super_dim}")
        self.name = args.name; print(f"Name: '{args.name}'")
        self.cmdargs = args; self.sampling_diagonal = (args.type == CFG_DIAG); self.theta = args.twist
        self.sampling_z = (args.type == CFG_Z); self.sampling_lc = (args.type == CFG_LC)
        self.cfg_grid_sz = None; self.calc_str = args.type; self.sampling = args.sampling
        self.no_ionic_step = args.noionicstep
        self.z_layer_sep = args.intspc
        if self.sampling_diagonal:
            print("Configuration sampling: DIAGONAL")
            self.cfg_grid_sz = DIAG_SAMPLE_LOW if args.sampling == 'low' else DIAG_SAMPLE_HIGH
        elif not self.sampling_z:
            print("Configuration sampling: GRID")
            self.cfg_grid_sz = GRID_SAMPLE_LOW if args.sampling == 'low' else GRID_SAMPLE_HIGH
        self.type_flag = TYPE_RELAX_BASIC
        if args.type in CFG_STRS:
            self.type_flag = TYPE_RELAX_CONFIG
        
        elif args.type == 'twist':
            self.type_flag = TYPE_TWISTED_CONFIG
            assert self.theta is not None, "Must provide a twist angle for twist calculations"
        elif args.type == 'norelax':
            self.type_flag = TYPE_NORELAX_BASIC
        self.run_relaxer = args.relax
        print(f"Relaxer: {'ON' if args.relax else 'OFF'}")
        assert (not args.relax) or (args.twist is not None), f"Must give a twist angle if running relaxer"
        
        print('Command line arguments initiated in InputData class constructor. Arguments: ', self.cmdargs)

        self.ROOT = args.dir
        self.do_vdW = args.vdw
        self.ediff0 = args.edinit if args.edinit is not None else EDIFF['relax']
        print(f"INITIAL EDIFF: {self.ediff0}")
        assert self.ediff0 is not None
        self.as_dfpt = args.dfpt
        self.fcut = args.fcut
        if self.fcut:
            print("USING: relaxation force cutoff (negative EDIFFG)")
        else:
            print("USING: relaxation energy cutoff (positive EDIFFG)")
        self.kpoints_is_gamma_centered = not args.mp
        self.calculation_list = tuple(dict.fromkeys(args.calc)) # Filter duplicates in calculation flag list
        self.input_imported = True
        self.__check_input_style()
        self.__parse_calculation_input()
        self.z = None; self.lcrange = None
        if self.sampling_z:
            print("Configuration sampling: INTERLAYER SPACING")
            assert os.path.isfile(self.ROOT + ISPC_SAMPLE_INAME), f"Path {self.ROOT + ISPC_SAMPLE_INAME} does not exist"
            with open(self.ROOT + ISPC_SAMPLE_INAME, 'r') as f:
                tup = tuple(map(float, f.read().splitlines()))
                assert len(tup) == 3 and tup[1] > tup[0] > 0 and tup[2] > 0, f"Invalid z input {tup}"
                self.z = np.linspace(tup[0], tup[1], num=int(tup[2]))
        elif self.sampling_lc:
            print("Configuration sampling: LATTICE CONSTANT")
            assert os.path.isfile(self.ROOT + LC_SAMPLE_INAME), f"Path {self.ROOT + LC_SAMPLE_INAME} does not exist"
            with open(self.ROOT + LC_SAMPLE_INAME, 'r') as f:
                tup = tuple(map(float, f.read().splitlines()))
                assert len(tup) == 3 and tup[1] > tup[0] > 0 and tup[2] > 0, f"Invalid LC input {tup}"
                self.lcrange = np.linspace(tup[0], tup[1], num=int(tup[2]))


        print('Final calculation list (except energies if specified):', self.calculation_list)

    # Confirm that the import is successful, exit otherwise
    def __check_import_status(self):
        if not self.input_imported:
            exit_with_error(ERR_NO_INPUT)

    # <Validators>
    def __check_type_flag(self):
        if self.type_flag not in TYPE_FLAGS:
            exit_with_error(ERR_BAD_TYPE_FLAG)
        # * 0: relax basic, 1: relax config, 2: no-relax basic (for phonon calculations)
        print('Calculation type flag validated.')

    def __check_root_dir(self):
        self.ROOT = checkPath(os.path.abspath(self.ROOT))
        if not os.path.isdir(self.ROOT):
            exit_with_error(ERR_BAD_INPUT_DIR)
        print('Root directory input validated.')
        
    def __check_vdW_settings(self):
        if self.do_vdW:
            print('Turning on van der Waals settings...')
        else:
            print("Using LDA electronic algorithm...")
        print('van der Waals flag validated.')
    
    def __check_kpoints_style(self):
        if self.kpoints_is_gamma_centered: # i.e. if user entered Gamma-centered kpoints
            print('Adopting Gamma-centered sampling of the Brillouin zone.')
        else:
            print('Adopting Monkhorst-Pack scheme of sampling the Brillouin zone.')
        print('Kpoints flag validated.')
    
    def __check_input_style(self):
        print('Validating command line inputs...')
        self.__check_import_status()
        self.__check_type_flag()
        self.__check_root_dir()
        self.__check_vdW_settings()
        self.__check_kpoints_style()
    
    # <Validators/>
    # Parse the input calculations and extract the energy calculation if there is one
    def __parse_calculation_input(self):
        print('Parsing list of calculations to perform...')
        self.__check_import_status()
        self.need_energies = False # whether or not to do total energy
        enerIndex = -1
        for i in self.calculation_list:
            if i not in CMD_LINE_ARG_LIST:
                exit_with_error(BAD_INPUT_ERR_MSG)
            if i == ENERGIES:
                self.need_energies = True
                enerIndex = self.calculation_list.index(ENERGIES)
        if PH in self.calculation_list and (PHDOS in self.calculation_list or PHBAND in self.calculation_list):
            exit_with_error(ERR_INCOMP_CALC1)
        if self.type_flag == TYPE_NORELAX_BASIC and (len(self.calculation_list) > 1 or ENERGIES not in self.calculation_list):
            exit_with_error(ERR_INCOMP_CALC2)
        if enerIndex >= 0:
            print('Recieved request to calculate total energies. Branching energy computation to different calculation pipeline...')
            self.calculation_list = list(self.calculation_list)
            self.calculation_list.pop(enerIndex)
            self.calculation_list = tuple(self.calculation_list)
            calcs = self.calculation_list if self.calculation_list else "(None)"
            print('Removed energy computation from main pipeline. Other calculations:', calcs)

    # Check the config type flag
    def get_type_flag(self):
        return self.type_flag
    
    # Query whether we are doing energy calculations or not
    def do_energy_calculation(self):
        return self.need_energies
    
    # Query whether there are any other calculations besides energy
    def do_nonenergy_calculations(self):
        if len(self.calculation_list) > 0:
            return True
        else:
            return False

    # Query whether we need to input a high-symmetry line kpoints object
    def need_line_kpoints(self):
        if (ELEBAND in self.calculation_list) or (PHBAND in self.calculation_list):
            return True
        else:
            return False

    def get_base_root_dir(self):
        return self.ROOT

    # Get calculation list
    def get_calculation_list(self):
        return list(self.calculation_list)

    # Get full calculation list, including energies
    def get_raw_calculation_list(self):
        if self.need_energies:
            return [ENERGIES] + self.get_calculation_list()
        else:
            return self.get_calculation_list()
        
    # Phonon supercell dimensions
    def get_super_dim(self):
        return self.super_dim

    def get_cfg_grid_sz(self):
        return self.cfg_grid_sz
    
    def sampling_is_diagonal(self):
        return self.sampling_diagonal

    def sampling_is_interlayer(self):
        return self.sampling_z
    
    def sampling_is_lc(self):
        return self.sampling_lc
    
    def get_tw_angle(self):
        return self.theta
    
    def get_interlayer_sampling(self):
        assert self.z is not None, "Calculation is not of interlayer sampling type"
        return self.z
    
    def get_lc_sampling(self):
        assert self.lcrange is not None, "Calculation is not of lattice constant sampling type"
        return self.lcrange

    def passname(self):
        return self.name

