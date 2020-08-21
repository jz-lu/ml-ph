from ___constants_names import PRGM_VERSION, PRGM_COOL_NAME
from datetime import datetime

def print_start_msg():
    print(PRGM_COOL_NAME)
    print('\n')
    print('Version:', PRGM_VERSION)
    print(datetime.now(), '\n')