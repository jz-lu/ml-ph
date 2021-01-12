from ____exit_with_error import exit_with_error

from ___constants_misc import *
from ___constants_names import CMD_LINE_ARG_LIST, ENERGIES, ELEBAND, PHBAND

from __directory_searchers import checkPath

import os

# Since this is a simple parsing class, we'll dispense with the usually reasonable move of mangling member variables to simulate privacy.
class InputData:
    input_imported = False

    def __init__(self, cmdline_arg_tuple): # cmdline_arg_tuple being an immutable tuple of sys.argv
        cmdline_arg_tuple = tuple(cmdline_arg_tuple)
        if (len(cmdline_arg_tuple) < 6) or (len(cmdline_arg_tuple) > 10):
            exit_with_error(GENERAL_ERR_USAGE_MSG) # Need at least one calculation in addition to settings
        self.cmdargs = cmdline_arg_tuple
        print('Command line arguments initiated in InputData class constructor. Arguments: ', self.cmdargs)
        try:
            self.type_flag = int(cmdline_arg_tuple[0])
        except ValueError as err:
            print('Error:', err)
            exit_with_error(ERR_BAD_TYPE_FLAG)
        self.__interlayer_distance = cmdline_arg_tuple[1]
        if self.__interlayer_distance <= 0:
            exit_with_error('Interlayer distance invalid.')

        self.ROOT = cmdline_arg_tuple[2]
        self.do_vdW = cmdline_arg_tuple[3]
        self.kpoints_is_gamma_centered = cmdline_arg_tuple[4]
        self.calculation_list = tuple(dict.fromkeys(cmdline_arg_tuple[5:])) # Filter duplicates in calculation flag list
        self.input_imported = True
        self.__check_input_style()
        self.__parse_calculation_input()
        if self.type_flag != 0:
            self.calculation_list = [ENERGIES] # In config space all we want are the energies, for now. TODO add dynamical matrix stuff?
        print('Final calculation list:', self.calculation_list)

    def __check_import_status(self):
        if not self.input_imported:
            exit_with_error(ERR_NO_INPUT)

    def __check_type_flag(self):
        if self.type_flag < 0:
            exit_with_error(ERR_BAD_TYPE_FLAG)
        print('Calculation type flag validated.')

    def __check_root_dir(self):
        self.ROOT = checkPath(self.ROOT)
        if not os.path.isdir(self.ROOT):
            exit_with_error(ERR_BAD_INPUT_DIR)
        print('Root directory input validated.')
        
    def __check_vdW_settings(self):
        if self.do_vdW == 'T' or self.do_vdW == 'True':
            self.do_vdW = True
            print('Turning on (SCAN-rVV10 functional, BPARAM=15.7) van der Waals settings...')
        elif self.do_vdW == 'F' or self.do_vdW == 'False':
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
        if enerIndex >= 0:
            print('Recieved request to calculate total energies. Branching energy computation to different calculation pipeline...')
            self.calculation_list = list(self.calculation_list)
            self.calculation_list.pop(enerIndex)
            self.calculation_list = tuple(self.calculation_list)
            print('Removed energy computation from main pipeline. New calculation list:', self.calculation_list)

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
        return self.calculation_list

    def get_interlayer_distance(self):
        return self.__interlayer_distance



