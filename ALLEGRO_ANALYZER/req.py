import os
import argparse

parser = argparse.ArgumentParser(description="DFT calculations on multilayered and twisted materials")
parser.add_argument("sampling", type=str, help="low or high", choices=['high', 'low'])
args = parser.parse_args()

for i in range(81 if args.sampling == 'high' else 9):
    pdir = f'shift_{i}/analyses/phonon'
    ndisp = len(list(filter(lambda x: x.startswith('disp'), os.listdir(pdir))))
    fails = []
    for j in range(1, ndisp+1):
        if len(os.listdir(f'shift_{i}/analyses/phonon/disp{j}')) <= 1:
            fails.append(j)
    if len(fails) > 0:
        os.chdir(pdir)
        fails = ",".join(fails)
        print(f"Resubmitting failed ph jobs: {fails}")
        stream = os.popen(f'sbatch --array={fails} EXECUTABLE_BAT_DNE')
        print(stream.read())
