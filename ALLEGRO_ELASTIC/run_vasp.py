# -*- coding: utf-8 -*-
"""
Created on Sun Nov 18 22:33:16 2018

@author: zoe
Modified by: jonathan

Submit jobs for wannierization
"""

import os 
import time 
import argparse


parser = argparse.ArgumentParser(description='Run VASP for elastic moduli calculations')
parser.add_argument('integers', metavar='N', type=int, nargs='+',
                    help='shear/bulk strain index')
parser.add_argument('folder', type=str, help='analysis folder (e.g. shear)')
parser.add_argument('-d', '--dir', type=str, help='working directory', default='.')
args = parser.parse_args()
ind = args.integers[0]
this_path = os.path.dirname(os.path.realpath(__file__))
print(f"This path is {this_path}")

folder = args.folder
main_dir = os.path.abspath(args.dir)
defdir = main_dir + '/' + folder + '/def' + str(ind) + '/'
os.chdir(defdir)
os.system('cp ' + main_dir + '/INCAR ' + ' ' + defdir)
os.system('cp ' + main_dir + '/POTCAR ' + ' ' + defdir)
os.system('cp ' + main_dir + '/KPOINTS ' + ' ' + defdir)
os.system('cp ' + this_path + '/vdw_kernel.bindat' + ' ' + defdir) # copy vasp vdw binary in

os.system('mpirun -np $SLURM_NTASKS /n/kaxiras_lab/vasp.5.4.4-centos7/vasp-O1.std')

if os.path.isfile('WAVECAR'):
    os.system('rm WAVECAR')

os.chdir(main_dir)

