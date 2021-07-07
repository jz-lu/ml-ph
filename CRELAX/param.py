import argparse
import os, sys
from gen_poscar_elastic import get_poscar_elastic, pathify
CODE_DIR = '/n/home04/jzlu/codes/ml-ph/CRELAX/'

folder_names = ['shear', 'stretch1', 'stretch2']

print(f"Args: {sys.argv}")
parser = argparse.ArgumentParser(description='Obtain energy calculations for elastic moduli fitting')
parser.add_argument("mat", type=str, help="material to analyze, e.g. MoS2, Gr, etc.")
parser.add_argument('-d', '--dir', type=str, help='working directory', default='.')
args = parser.parse_args()
ROOT = pathify(args.dir)
for i in ['INCAR', 'POTCAR', 'KPOINTS']:
    assert os.path.isfile(args.dir + '/' + i), f"Missing file {i} in {args.dir}"
print(f"Args: {sys.argv[1:]}")

# Run the file maker
n_eta, updated_root = get_poscar_elastic(args)
ROOT = pathify(updated_root)

# Run VASP
for folder in folder_names:
    fpath = ROOT + folder + '/' + 'ELASTIC_BAT_DNE'
    with open(fpath, 'w') as f:
        f.write('#!/bin/bash\n')
        f.write(f'#SBATCH -J {folder}\n')
        f.write('#SBATCH -n 16\n#SBATCH -t 8:00:00\n')
        f.write('#SBATCH -p kaxiras,shared\n#SBATCH --mem-per-cpu=4000\n')
        f.write(f'#SBATCH -o {folder}_%a_%A.out\n#SBATCH -e er_{folder}_%a_%A.err')
        f.write('#SBATCH --mail-type=END,FAIL\n#SBATCH --mail-user=jlu@college.harvard.edu\n')
        f.write('module load intel/17.0.4-fasrc01\n')
        f.write('module load impi/2017.2.174-fasrc01\n')
        f.write('module load python\n')
        f.write('module list\n')
        f.write('export TASKID=$SLURM_ARRAY_TASK_ID\n')
        f.write('echo "Task ID: " $TASKID\n')
        f.write(f'python {CODE_DIR}/run_vasp.py $TASKID\n')
    cmd = f'sbatch --array=0-{n_eta-1} {fpath}'
    print(f"Running command `{cmd}` to shell...")
    stream = os.popen(cmd)
    print(stream.read())
print("All jobs submitted. Run analyzer (`read_elastic.py`) upon completion.")



