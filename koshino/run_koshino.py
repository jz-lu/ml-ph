"""
Streamlined Koshino model calculations over many twist angles.
"""
import os, argparse
import numpy as np
RELATIVE_PATH = 'ml-ph/koshino/phonon_2D.jl'
LOCAL_PATH = '/Users/jonathanlu/Documents/'
CLUSTER_PATH = '/n/home04/jzlu/codes/'

def get_G_Cutoff(theta):
    if theta >= 6:
        return 7
    if 3 <= theta < 6:
        return 9
    if 1 <= theta < 3:
        return 11
    else:
        return 13

parser = argparse.ArgumentParser(description="Koshino model calculations over many twist angles")
parser.add_argument("-N", "--gridsz", type=int, help='enter N for NxN sampling mesh', default=51)
parser.add_argument("-o", "--out", type=str, help='output path', default='.')
parser.add_argument("-m", "--mesh", action='store_true', help="run DOS over mesh instead of band")
parser.add_argument("-c", "--cluster", action='store_true', help="run on cluster instead of local")
parser.add_argument("range", nargs=3, help="start, end, step (in degrees)")
args = parser.parse_args()

rng = list(map(float, args.range))
rng[1] += rng[2]
thetas = np.arange(*rng)
print(f"Twist angles: {thetas}")
do_DOS = '--mesh' if args.mesh else ''
path = (CLUSTER_PATH if args.cluster else LOCAL_PATH) + RELATIVE_PATH

for theta in thetas:
    theta = round(theta, 4)
    print(f"WORKING ON: {theta} deg...", flush=True)
    cmd_here = f"julia {path} -d {theta} -N {args.gridsz} -o {args.out} -c {get_G_Cutoff(theta)} {do_DOS}"
    print(f"RUNNING: `{cmd_here}`", flush=True)
    stream = os.popen(cmd_here)
    print(stream.read())
    print("DONE", flush=True)
print("All done!")