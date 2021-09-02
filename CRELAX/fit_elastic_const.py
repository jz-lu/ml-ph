"""
Script for analyzing output of elastic moduli calculations kickedd off in `param.py`.
Python equivalent of `fit_elastic_const.m`
"""
import numpy as np
import argparse
import os

parser = argparse.ArgumentParser(description="Output analyzer for elastic moduli calculations")
parser.add_argument("-n", "--name", type=str, help='calculation name', default='')
parser.add_argument("-t", "--type", type=str, help="basic, config, diag, norelax, z, lc, or twist",
                    default='basic', choices=TYPE_STRS)
parser.add_argument("-s", "--sampling", type=str, help='low or high', default='low', choices=['low', 'high'])
parser.add_argument("-e", "--edinit", type=float, help="initial EDIFF for TLT algorithm", default=None)
parser.add_argument("--twist", type=float, help="give a twist angle", default=None)
parser.add_argument("-d", "--dir", type=str, help="main working directory", default='.')
parser.add_argument("-v", "--vdw", action="store_true", help="use van der Waals corrections")
parser.add_argument("-f", "--fcut", action="store_true", help="use force cutoff instead of energy for EDIFFG")
parser.add_argument("--dfpt", action="store_true", help="use DFPT instead of frozen phonon")
parser.add_argument("-m", "--mp", action="store_true", help="use for MP k-points mesh, default: Gamma")
parser.add_argument("-r", "--relax", action="store_true", help="use relaxed (non-uniform) configurations")
parser.add_argument("--super", type=int, help="phonon supercell", default=SUPER_DIM[0])
parser.add_argument("calc", nargs="+", help="calculations list: energies, ph, eleband, eledos", default=[ENERGIES])
cmdargs = parser.parse_args()
