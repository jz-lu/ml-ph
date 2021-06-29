from ____exit_with_error import exit_with_error
from ___constants_vasp import BPARAM
from ___constants_misc import *
from ___constants_config import GRID_SAMPLE_HIGH, GRID_SAMPLE_LOW
from ___constants_names import (
    CMD_LINE_ARG_LIST, 
    ENERGIES, ELEBAND, PHBAND, PHDOS, PH, 
    TYPE_FLAGS, 
    TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC
)
from __directory_searchers import checkPath
import os

# Since this is a simple parsing class, we'll dispense with the usually reasonable move of mangling member variables to simulate privacy.
class InputData:
    input_imported = False

    def __init__(self, cmdline_arg_tuple):
        cmdline_arg_tuple = tuple(cmdline_arg_tuple)
        if (len(cmdline_arg_tuple) < 5) or (len(cmdline_arg_tuple) > 10):
            exit_with_error(GENERAL_ERR_USAGE_MSG) # Need at least one calculation in addition to settings
        self.cmdargs = cmdline_arg_tuple
        print('Command line arguments initiated in InputData class constructor. Arguments: ', self.cmdargs)
        try:
            self.cfg_grid_sz = None
            if cmdline_arg_tuple[0] in ['LOW', 'HIGH', 'L', 'H']:
                self.cfg_grid_sz = GRID_SAMPLE_LOW if cmdline_arg_tuple[0] in ['LOW', 'L'] else GRID_SAMPLE_HIGH
                cmdline_arg_tuple[0] = 1
            self.type_flag = int(cmdline_arg_tuple[0])
        except ValueError as err:
            print('Error:', err)
            exit_with_error(ERR_BAD_TYPE_FLAG)

        self.ROOT = cmdline_arg_tuple[1]
        self.do_vdW = cmdline_arg_tuple[2]
        self.kpoints_is_gamma_centered = cmdline_arg_tuple[3]
        self.calculation_list = tuple(dict.fromkeys(cmdline_arg_tuple[4:])) # Filter duplicates in calculation flag list
        self.input_imported = True
        self.__check_input_style()
        if len(self.calculation_list) == 0:
            self.calculation_list = [ENERGIES]
            # TODO add dynamical matrix stuff?
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
        if self.do_vdW == 'T' or self.do_vdW == 'True':
            self.do_vdW = True
            print("Using SCAN electronic algorithm...")
            print(f'Turning on (SCAN-rVV10 functional, BPARAM={BPARAM}) van der Waals settings...')
        elif self.do_vdW == 'F' or self.do_vdW == 'False':
            print("Using LDA electronic algorithm...")
            self.do_vdW = False
        else:
            print(ERR_INVALID_VDW_FLAG)
            exit_with_error(GENERAL_ERR_USAGE_MSG)
        print('van der Waals flag validated.')
    
    def __check_kpoints_style(self):
        if (self.kpoints_is_gamma_centered).lower()[0] == 'g': # i.e. if user entered Gamma-centered kpoints
            self.kpoints_is_gamma_centered = True
            print('Adopting Gamma-centered sampling of the Brillouin zone.')
        else:
            self.kpoints_is_gamma_centered = False
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

