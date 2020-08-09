import sys
from datetime import datetime

def exit_with_error(errorMsg):
    print(errorMsg)
    sys.exit('Time program terminated:', datetime.now())