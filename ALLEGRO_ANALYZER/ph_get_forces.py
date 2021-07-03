# Extract FORCE_SETS from VASP calculations using phonopy
import os, sys, copy
from ___helpers_parsing import update, succ, warn, err, is_flag, check_not_flag
from ___constants_names import (
    PHDISP_STATIC_NAME, PH_FORCE_SETS_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, CONFIG_SUBDIR_NAME
)
from __directory_searchers import findDirsinDir, checkPath
from __dirModifications import build_dir
from __ph_processing import ph_generate_forcesets
from __dirModifications import move

def get_many_forces(indir, twist=False):
    indir = checkPath(indir)
    assert os.path.isdir(indir), f"Directory {indir} does not exist"
    layers = []; configs = None
    if twist:
        layers = findDirsinDir(indir, MONOLAYER_DIR_NAME, searchType='start')
        layers = [build_dir([indir, layer, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME]) for layer in layers]
        configs = findDirsinDir(indir + CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME, searchType='start')
        configs = [build_dir([indir, CONFIG_DIR_NAME, config, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME]) for config in configs]
    else:
        configs = findDirsinDir(indir, CONFIG_SUBDIR_NAME, searchType='start')
        configs = [build_dir([indir, config, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME]) for config in configs]
    paths = layers + configs # concatenate all paths together
    for path in paths: # intralayer terms
        print(f"Generating {PH_FORCE_SETS_NAME} in {path}...")
        assert os.path.isdir(path), f"Directory {path} does not exist"
        disps = findDirsinDir(path, PHDISP_STATIC_NAME, searchType='start')
        ph_generate_forcesets(path, len(disps), path_pad=ANALYSIS_DIR_NAME)
        if not os.path.isfile(path + PH_FORCE_SETS_NAME):
            err(f"Error: could not find {PH_FORCE_SETS_NAME} in {path}. Check phonopy output for log.")
    return

if __name__ == '__main__':
    args = copy.deepcopy(sys.argv)[1:]; i = 0; n = len(args)
    indir = '.'; outdir = None; ftype = 'basic'
    USAGE_ERR_MSG = f"Usage: python3 {sys.argv[0]} -dir <MAIN PHONON DIRECTORY> -o <PRINT DIR (optional)> --tw (if twisted) --cfg (if config)"
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-dir':
            i += 1; check_not_flag(args[i]); indir = checkPath(args[i]); i += 1
        elif args[i] == '-o':
            i += 1; check_not_flag(args[i]); outdir = checkPath(args[i]); i += 1
        elif args[i] == '--tw':
            ftype = 'twist'; i += 1
        elif args[i] == '--cfg':
            ftype = 'config'; i += 1
        elif args[i] == '--usage':
            print(USAGE_ERR_MSG)
            sys.exit(0)
        else:
            err(USAGE_ERR_MSG)
    outdir = indir if outdir is None else outdir
    indir = checkPath(os.path.abspath(indir)); outdir = checkPath(os.path.abspath(outdir))

    if ftype == 'twist':
        update("Mode: twisted")
        get_many_forces(indir, twist=True)
        succ(f"Successfully generated intra- (and inter-)layer force constants files in {indir}")
    elif ftype == 'config':
        update("Mode: configuration")
        get_many_forces(indir, twist=False)
        succ(f"Successfully generated per-configuration force constants files in {indir}")
    else:
        # Call FORCE_SETS generator
        update("Mode: basic")
        disps = findDirsinDir(indir, PHDISP_STATIC_NAME, searchType='start')
        ph_generate_forcesets(indir, len(disps), path_pad=ANALYSIS_DIR_NAME)
        if not os.path.isfile(indir + PH_FORCE_SETS_NAME):
            err(f"Error: could not find {indir + PH_FORCE_SETS_NAME}. Check phonopy output for log.")
        if indir != outdir:
            # Move the FORCE_SETS file to the desired directory
            move(PH_FORCE_SETS_NAME, indir, newPath=outdir)
        succ("Successfully generated force constants file")



