"""
A small script to make some POSCAR adjustments. The most interesting
one is rotating atoms in a unit cell.
"""

from pymatgen.io.vasp.inputs import Poscar
from pymatgen.core.structure import Structure
import numpy as np
import numpy.linalg as LA
import argparse, os
from ___constants_names import POSCAR_NAME
from __directory_searchers import checkPath
from ___helpers_parsing import succ
from math import cos, sin

def rotate_atoms(theta, p, out=None):
    # everything is row-wise, so a little algebra shows how to minimize transposes
    R = np.array([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0,0,1]])
    pos = p.structure.cart_coords @ R.T
    A0 = p.structure.lattice.matrix
    species = [s.name for s in p.structure.species]
    struct = Structure(A0, species, pos, coords_are_cartesian=True)
    if out is not None:
        struct.to(fmt='poscar', filename=out) # write unstreched POSCAR to file 
        succ(f"Wrote updated POSCAR with rotation {np.rad2deg(theta)} to {out}")
    return struct

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Adjust POSCARs, such as rotating atoms")
    parser.add_argument("deg", type=float, help="rotation angle", default=None)
    parser.add_argument("-d", "--dir", type=str, help="directory containing POSCAR", default='.')
    parser.add_argument("-p", "--pos", type=str, help="POSCAR name", default=POSCAR_NAME)
    parser.add_argument("-o", "--out", type=str, help="Output POSCAR name", default=POSCAR_NAME+"_out")
    args = parser.parse_args()

    args.dir = checkPath(os.path.abspath(args.dir))
    ppath = args.dir + args.pos
    assert os.path.isfile(ppath), f"Invalid path {ppath}"
    p = Poscar.from_file(ppath)
    rotate_atoms(np.deg2rad(args.deg), p, out=args.dir+args.out)

