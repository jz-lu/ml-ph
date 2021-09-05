"""
Script for removing all files in a twisted calculations 
not necessary to construct dynamical matrix.
"""
import os
import argparse
from __directory_searchers import checkPath
from ___helpers_parsing import update, succ

parser = argparse.ArgumentParser(description="Remove unnecessary files for twisted phonon calcs")
parser.add_argument("-d", "--dir", type=str, help="root calc directory", default='.')
args = parser.parse_args()
ROOT = checkPath(os.path.abspath(args.dir))
PREFIXES = tuple(['FORCE', 'POSCAR', 'CONTCAR'])
SUFFIXES = tuple(['.txt'])

update("Cleaning files...")
for root, dirs, files in os.walk(ROOT, topdown=True):
   for name in files:
       trash = False
       if name.startswith(PREFIXES) or name.endswith(SUFFIXES):
            path = os.path.join(root, name)
            print(f"REMOVING: {path}")
            os.remove(path)

update("Clearning folders...")
for root, dirs, files in os.walk(ROOT, topdown=True):
    for name in dirs:
        path = os.path.join(root, name)
        if not os.listdir(path):
            print(f"CLEARING: {path}")
            os.rmdir(path)

succ("Cleaning complete.")
   
