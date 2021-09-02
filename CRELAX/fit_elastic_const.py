"""
Script for analyzing output of elastic moduli calculations kickedd off in `param.py`.
Python equivalent of `fit_elastic_const.m`
"""
import numpy as np
import argparse
import os
from param import folder_names

parser = argparse.ArgumentParser(description="Output analyzer for elastic moduli calculations")
parser.add_argument("-d", "--dir", type=str, help="main working directory", default='.')
parser.add_argument("lconst", type=str, help="lattice constant")
args = parser.parse_args()
a0 = args.lconst

os.chdir(args.dir)
infiles = [name + '.txt' for name in folder_names]
assert set(infiles).issubset(set(os.listdir())), f"Input files not found, run `read_elastic.py` first!"

