import os
import numpy as np
import argparse
from consts import CODE_DIR

parser = argparse.ArgumentParser(description="Make band plots")
parser.add_argument("-d", "--dir", type=str, help='main directory path', default=".")
parser.add_argument("-m", "--data", type=str, help='data directory path', default="..")
parser.add_argument("-n", "--name", type=str, help='material name', default="(mat)")
parser.add_argument("-r", "--relax", action="store_true", help='run relaxer')
parser.add_argument("range", nargs=3, type=int, help="theta: (start, end, number of)")
args = parser.parse_args()

t_start, t_end, n_t = tuple(args.range)
thetas = np.linspace(t_start, t_end, n_t)
main_dir = os.path.abspath(args.dir)
data_dir = os.path.abspath(args.data)
r = "-r" if args.relax else ""
if r == "-r":
    print("RUNNING: relaxer", flush=True)
os.chdir(main_dir)

for theta in thetas:
    theta = round(theta, 6)
    print("RUNNING: theta =", theta, flush=True)
    os.mkdir(f't{theta}')
    print(os.popen(f"python {CODE_DIR} -d {data_dir} -o t{theta} -n {args.name} {r} {theta}").read())