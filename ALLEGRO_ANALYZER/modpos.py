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

# Change the lattice constant
def update_lc(p, lc, out=None):
    p = p.as_dict()
    old_mat = np.array(p['structure']['lattice']['matrix'])
    old_lc = LA.norm(old_mat[0])
    m = old_mat[-1,-1] * np.eye(3)
    m[:2,:2] = lc/old_lc * old_mat[:2,:2]
    p['structure']['lattice']['matrix'] = m
    p = Poscar.from_dict(p)
    if out is not None:
        p.write_file(out)
        succ(f"Wrote updated POSCAR with lattice constant {lc} to {out}")
    return p

# Change the vacuum spacing
def update_vacspc(p, z, out=None):
    n_at = len(p.structure.species)
    old_z = p.structure.lattice.matrix[-1,-1]
    p = p.as_dict()
    for i in range(n_at):
        p['structure']['sites'][i]['abc'][-1] *= z/old_z
    p['structure']['lattice']['matrix'][-1][-1] = z
    p = Poscar.from_dict(p)
    if out is not None:
        p.write_file(out)
        succ(f"Wrote updated POSCAR with vacuum spacing {z} to {out}")
    return p

def rotate_atoms(theta, p, out=None): # theta in rad
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
    parser = argparse.ArgumentParser(description="Adjust POSCARs, such as rotating atoms or changing constants")
    parser.add_argument("type", choices=['lc', 'z', 'rot'])
    parser.add_argument("param", type=float, help="parameter (e.g. angle for rotation, lattie constant, etc.", default=None)
    parser.add_argument("-d", "--dir", type=str, help="directory containing POSCAR", default='.')
    parser.add_argument("-p", "--pos", type=str, help="POSCAR name", default=POSCAR_NAME)
    parser.add_argument("-o", "--out", type=str, help="Output POSCAR name", default=POSCAR_NAME+"_out")
    args = parser.parse_args()

    args.dir = checkPath(os.path.abspath(args.dir))
    ppath = args.dir + args.pos; outpath = args.dir+args.out
    assert os.path.isfile(ppath), f"Invalid path {ppath}"
    p = Poscar.from_file(ppath)
    if args.type == 'rotate':
        rotate_atoms(np.deg2rad(args.param), p, out=outpath)
    elif args.type == 'lc':
        update_lc(p, args.param, out=outpath)
    elif args.type == 'z':
        update_vacspc(p, args.param, out=outpath)


