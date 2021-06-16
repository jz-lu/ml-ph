# Collection of useful functions to help parse commandline input
import sys

class bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'

def succ(s):
    print(bcolors.OKGREEN + s + bcolors.ENDC)

def warn(s):
    print(bcolors.WARNING + s + bcolors.ENDC)

def err(s):
    print(bcolors.FAIL + s + bcolors.ENDC)
    sys.exit(1)

def is_flag(s):
    return s[0] == '-'

def check_not_flag(s):
    if is_flag(s):
        err(f'Error: flag "{s}" found before argument for previous flag given')