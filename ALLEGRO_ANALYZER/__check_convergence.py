# Checks whether the relaxation has converged or not
from ___constants_vasp import VASP_OUTFILE_CONVERGENCE_LINE_CUT, NSW
from ___constants_names import OSZICAR_NAME
from __directory_searchers import checkPath

# Function to read last N lines of the file  
def last_N_lines(fname, N): 
    # opening file using with() method so cleanup is automatic
    outstr = ''
    with open(fname) as f: 
        for line in (f.readlines() [-N:]): 
            outstr += line
    return outstr

# Uses last_N_lines to check convergence of file
def check_if_not_converged(dirName, fileName):
    print('Now checking for convergence...')
    dirName = checkPath(dirName); filePath = dirName + fileName

    # Flag keywords we'll look for
    convergence_err_keywords = ['fatal error in bracketing', 'copy CONTCAR']

    # Check the file for convergence errors
    string = last_N_lines(filePath, VASP_OUTFILE_CONVERGENCE_LINE_CUT)

    fatal_bracket = True
    for i in convergence_err_keywords:
        if string.find(i) == -1:
            fatal_bracket = False
            break
    nsw = None; nsw_reached = False
    with open(dirName + OSZICAR_NAME, 'r') as f:
        nsw = f.readlines()[-1].split(' ')
        while not nsw[0]:
            del nsw[0]
        nsw = int(nsw[0])
        print(f"Took {nsw} out of max {NSW['relax']} steps")
        nsw_reached = True if nsw >= NSW['relax'] else False
    
    return fatal_bracket or nsw_reached
