# Checks whether the relaxation has converged or not
from ___constants_vasp import VASP_OUTFILE_CONVERGENCE_LINE_CUT

from __directory_searchers import checkPath

# Function to read last N lines of the file  
def last_N_lines(fname, N): 
    # opening file using with() method so cleanup is automatic
    outstr = ''
    with open(fname) as file: 
        for line in (file.readlines() [-N:]): 
            outstr += line
            #print(line, end ='')
    
    return outstr

# Uses last_N_lines to check convergence of file
def check_if_not_converged(dirName, fileName):
    print('Now checking for convergence...')
    dirName = checkPath(dirName)
    filePath = dirName + fileName

    # Flag keywords we'll look for
    convergence_err_keywords = ['fatal error in bracketing', 'copy CONTCAR']

    # Check the file for convergence errors
    string = last_N_lines(filePath, VASP_OUTFILE_CONVERGENCE_LINE_CUT)

    not_converged_yet = True # Initialized to true, if we don't find one of those keywords we falsify it
    for i in convergence_err_keywords:
        if string.find(i) == -1:
            not_converged_yet = False
            break
    
    return not_converged_yet
