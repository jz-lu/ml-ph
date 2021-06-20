# Extract FORCE_SETS from VASP calculations using phonopy
import os, sys, copy
from ___helpers_parsing import succ, warn, err, is_flag, check_not_flag
from ___constants_names import PHDISP_STATIC_NAME, PH_FORCE_SETS_NAME
from __directory_searchers import findDirsinDir, checkPath
from __ph_processing import ph_generate_forcesets
from __dirModifications import move

if __name__ == '__main__':
    # Parse input
    args = copy.deepcopy(sys.argv)[1:]; i = 0; n = len(args)
    indir = '.'
    outdir = '.'
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-dir':
            i += 1; check_not_flag(args[i]); indir = checkPath(args[i]); i +=1
        if args[i] == '-o':
            i += 1; check_not_flag(args[i]); outdir = checkPath(args[i]); i +=1
        else:
            err(f"Usage: python3 {sys.argv[0]} -dir <DISP DIRECTORY> -o <PRINT DIR (optional)>")

    # Call FORCE_SETS generator
    disps = findDirsinDir(indir, PHDISP_STATIC_NAME, searchType='start')
    ndisp = '%03d'%(len(disps))
    ph_generate_forcesets(indir, ndisp)
    if PH_FORCE_SETS_NAME not in os.listdir(indir):
        err(f"Error: force constants file not found at {indir}. Check phonopy output for log.")
    if indir != outdir:
        # Move the FORCE_SETS file to the desired directory
        move(PH_FORCE_SETS_NAME, indir, newPath=outdir)
    succ("Successfully generated force constants file")
else:
    err(f"File {__name__} has no function when imported")

