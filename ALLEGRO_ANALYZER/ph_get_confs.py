from pymatgen.io.vasp.inputs import Poscar, Kpoints
from __get_ph_analysis import ph_create_band_conf, ph_create_mesh_conf
from __directory_searchers import checkPath
from ___constants_names import KPOINTS_NAME, KPOINTS_LINE_NAME, POSCAR_NAME
from ___helpers_parsing import greet, succ, err, warn, is_flag, check_not_flag
import os, sys

if __name__ == '__main__':
    kdir = '.'; pdir = '.'; outdir = None; conf_t = 'm'; args = sys.argv[1:]; i = 0; n = len(args)
    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-pdir':
            i += 1; check_not_flag(args[i]); pdir = checkPath(args[i]); i += 1
        elif args[i] == '-kdir':
            i += 1; check_not_flag(args[i]); kdir = checkPath(args[i]); i += 1
        elif args[i] == '-o':
            i += 1; check_not_flag(args[i]); outdir = checkPath(args[i]); i += 1
        elif args[i] == '-ctype':
            i += 1; check_not_flag(args[i]); assert args[i] in ['b', 'm'], "-ctype flag must be 'b' or 'm'"
            conf_t = args[i]; i += 1
        else:
            err(f"Usage: python3 {sys.argv[0]} -dir <DISP DIRECTORY> -o <PRINT DIR (optional)>")
    if outdir is None:
        outdir = pdir
    kname = KPOINTS_LINE_NAME if conf_t == 'b' else KPOINTS_NAME
    assert os.path.isfile(kdir + kname), f"{kdir} does not contain {kname}"
    assert os.path.isfile(pdir + POSCAR_NAME), f"{pdir} does not contain {POSCAR_NAME}"
    p = Poscar.from_file(pdir + POSCAR_NAME); k = Kpoints.from_file(kdir + kname)
    greet(f"Building phonopy {'band' if conf_t == 'b' else 'mesh'} configuration file...")

    if conf_t == 'b':
        ph_create_band_conf(k, p, outdir)
    else:
        ph_create_mesh_conf(k, p, outdir)
    succ(f"Successfully wrote configuration file to {outdir}")

else:
    err(f"Error: program {sys.argv[0]} cannot be imported")