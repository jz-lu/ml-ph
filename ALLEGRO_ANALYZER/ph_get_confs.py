from pymatgen.io.vasp.inputs import Poscar, Kpoints
from __get_ph_analysis import ph_create_band_conf, ph_create_mesh_conf
from __directory_searchers import checkPath
from ___constants_names import KPOINTS_NAME, KPOINTS_LINE_NAME, POSCAR_NAME
from ___helpers_parsing import greet, succ, err, warn, is_flag, check_not_flag
import os, sys

if __name__ == '__main__':
    kpath = '.'; ppath = '.'; outdir = None; conf_t = 'm'; args = sys.argv[1:]; i = 0; n = len(args)
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-pp':
            i += 1; check_not_flag(args[i]); ppath = args[i]; i += 1
        elif args[i] == '-kp':
            i += 1; check_not_flag(args[i]); kpath = args[i]; i += 1
        elif args[i] == '-o':
            i += 1; check_not_flag(args[i]); outdir = checkPath(args[i]); i += 1
        elif args[i] == '-ctype':
            i += 1; check_not_flag(args[i]); assert args[i] in ['b', 'm'], "-ctype flag must be 'b' or 'm'"
            conf_t = args[i]; i += 1
        else:
            err(f"Usage: python3 {sys.argv[0]} -ctype <m or b> -pdir <relaxed POSCAR DIR> -kdir <KPOINTS dir> -o <PRINT DIR>")
    if outdir is None:
        outdir = ppath
    kname = KPOINTS_LINE_NAME if conf_t == 'b' else KPOINTS_NAME
    if os.path.isdir(ppath):
        ppath = checkPath(ppath) + POSCAR_NAME
    if os.path.isdir(kpath):
        kpath = checkPath(kpath) + kname
    assert os.path.isfile(kpath), f"{kpath} does not contain {kname}"
    assert os.path.isfile(ppath), f"{ppath} does not contain {POSCAR_NAME}"
    p = Poscar.from_file(ppath); k = Kpoints.from_file(kpath)
    greet(f"Building phonopy {'band' if conf_t == 'b' else 'mesh'} configuration file...")

    if conf_t == 'b':
        ph_create_band_conf(k, p, outdir)
    else:
        ph_create_mesh_conf(k, p, outdir)
    succ(f"Successfully wrote configuration file to {outdir}")

else:
    err(f"Error: program {sys.argv[0]} cannot be imported")