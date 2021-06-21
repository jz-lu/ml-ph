import sys
from datetime import datetime
from ___helpers_parsing import err

def exit_with_error(errorMsg):
    err(errorMsg, q='Time program terminated: ' + str(datetime.now()))
