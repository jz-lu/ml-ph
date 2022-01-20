"""
Script for removing all files in a twisted calculation;
not necessary to construct dynamical matrix.
"""
import os
import argparse
from __directory_searchers import checkPath
from ___helpers_parsing import update, succ

parser = argparse.ArgumentParser(description="Remove unnecessary files for twisted phonon calcs")
parser.add_argument("-d", "--dir", type=str, help="root calc directory", default='.')
parser.add_argument("-v", "--verbose", action="store_true", help="output a log of all removals")
args = parser.parse_args()
ROOT = checkPath(os.path.abspath(args.dir))
PREFIXES = tuple(['FORCE', 'SPOSCAR', 'POSCAR', 'CONTCAR', 'INCAR'])
SUFFIXES = tuple(['.txt'])

def vprint(s):
    if args.verbose:
        print(s, flush=True)

update("Cleaning files...")
for root, dirs, files in os.walk(ROOT, topdown=True):
   for name in files:
       if not (name.startswith(PREFIXES) or name.endswith(SUFFIXES)):
            path = os.path.join(root, name)
            vprint(f"REMOVING: {path}")
            os.remove(path)

update("Cleaning folders...")
for root, dirs, files in os.walk(ROOT, topdown=True):
    for name in dirs:
        path = os.path.join(root, name)
        if not os.listdir(path):
            vprint(f"CLEARING: {path}")
            os.rmdir(path)

succ("Cleaning complete.")
   
