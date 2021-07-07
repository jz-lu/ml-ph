from ____exit_with_error import exit_with_error
from ___constants_vasp import BPARAM
from ___constants_misc import *
from ___constants_config import GRID_SAMPLE_HIGH, GRID_SAMPLE_LOW, DIAG_SAMPLE_HIGH, DIAG_SAMPLE_LOW
from ___constants_names import (
    CMD_LINE_ARG_LIST, 
    ENERGIES, ELEBAND, PHBAND, PHDOS, PH, 
    TYPE_FLAGS, 
    TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC
)
from __directory_searchers import checkPath
import os, argparse

# Since this is a simple parsing class, we'll dispense with the usually reasonable move of mangling member variables to simulate privacy.
class InputData:
    input_imported = False

    def __init__(self, args):
        self.cmdargs = args; self.sampling_diagonal = (args.type == 'diag'); self.theta = args.twist
        self.cfg_grid_sz = None; self.calc_str = args.type; self.sampling = args.sampling
        if self.sampling_diagonal:
            self.cfg_grid_sz = DIAG_SAMPLE_LOW if args.sampling == 'low' else DIAG_SAMPLE_HIGH
        else:
            self.cfg_grid_sz = GRID_SAMPLE_LOW if args.sampling == 'low' else GRID_SAMPLE_HIGH
        self.type_flag = TYPE_RELAX_BASIC
        if args.type in ['config', 'diag']:
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
        self.kpoints_is_gamma_centered = not args.mp
        self.calculation_list = tuple(dict.fromkeys(args.calc)) # Filter duplicates in calculation flag list
        self.input_imported = True
        self.__check_input_style()
        self.__parse_calculation_input()

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

    def get_cfg_grid_sz(self):
        return self.cfg_grid_sz
    
    def sampling_is_diagonal(self):
        return self.sampling_diagonal
    
    def get_tw_angle(self):
        return self.theta

