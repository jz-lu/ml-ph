import argparse
import numpy as np
import os
from consts import ALLEGRO_DIR

parser = argparse.ArgumentParser(description="Submit band plots job to cluster")
parser.add_argument("-d", "--dir", type=str, help='main directory path', default=".")
parser.add_argument("-m", "--data", type=str, help='data directory path', default="..")
parser.add_argument("-n", "--name", type=str, help='material name', default="(mat)")
parser.add_argument("-r", "--relax", action="store_true", help='run relaxer')
parser.add_argument("range", nargs=3, type=int, help="theta: (start, end, number of)")
args = parser.parse_args()

t_start, t_end, n_t = tuple(args.range)
thetas = np.linspace(t_start, t_end, n_t)
main_dir = os.path.abspath(args.dir)
path = main_dir + '/' + 'batband'
r = "-r" if args.relax else ""

with open(path, 'w') as f:
    f.write("#!/bin/bash\n")
    f.write(f"#SBATCH --job-name=rs-{args.name}\n")
    f.write("#SBATCH -N 1\n#SBATCH -n 24\n")
    f.write("#SBATCH -t 16:00:00\n")
    f.write("#SBATCH -p kaxiras,shared\n")
    f.write("#SBATCH --mem-per-cpu=4000\n")
    f.write("#SBATCH -o rs_%a.out\n#SBATCH -e er_rs_%a.err\n")
    f.write("#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=jlu@college.harvard.edu\n")
    f.write("source activate $HOME/anaconda_env\n")
    f.write(f'WDIR="{main_dir}' + '"\necho "WD: ${WDIR}"\n')
    f.write(f'ALLEGRO_DIR="{ALLEGRO_DIR}"\n')
    f.write("module load julia\nmodule list\nsource activate $HOME/anaconda_env\n")
    f.write(f'echo "Starting calculations..."\npython3 make_bands.py -d {args.dir} --data {args.data} {r} -n {args.name} {t_start} {t_end} {n_t}\necho "Calculations complete!"')

print(f"Submitting {path}...", flush=True)
print(os.popen(f"sbatch {path}").read(), flush=True)















