import os, shutil
import argparse

parser = argparse.ArgumentParser(description="Requeue failed jobs on cluster")
parser.add_argument("sampling", type=str, help="low or high", choices=['high', 'low'])
parser.add_argument("-s", "--skip", nargs="+", help="numbers to skip", default=[])
parser.add_argument("-c", "--cfg", action="store_true", help="configuration requeue")
args = parser.parse_args()

d = os.getcwd()
skip = list(map(lambda x: int(x), args.skip))
print(f"Skipping: {skip}")
N = 81 if args.sampling == 'high' else 9
if args.cfg:
    fails = []
    for i in range(N):
        if i in skip:
            continue
        s1 = os.popen(f"grep 'not available' shift_{i}/relaxation/relax.err | wc -l").read()
        s2 = os.popen(f"grep 'rror' shift_{i}/relaxation/relax.err | wc -l").read()
        if int(s1) > 0 or int(s2) > 0:
            shutil.rmtree(f'shift_{i}/relaxation/')
            shutil.rmtree(f'shift_{i}/analyses/')
            fails.append(i)
    fails = ",".join(list(map(lambda x: str(x), fails)))
    print(f"Requeueing fails: {fails}")
    stream = os.popen(f"sbatch --array={fails} EXECUTABLE_BAT_DNE")
    print(stream.read())
    
else:
    no_fails = True
    for i in range(N):
        if i in skip:
            continue
        pdir = f'shift_{i}/analyses/phonon/'
        ndisp = len(list(filter(lambda x: x.startswith('disp'), os.listdir(pdir))))
        fails = []
        for j in range(1, ndisp+1):
            subdir = pdir + f'disp{j}/'
            never_started = len(os.listdir(subdir)) <= 1
            if never_started:
                print(f"[{i}] Found job never started at disp{j}")
                fails.append(j)
                continue
            subdir += 'analyses/'
            bad_term = int(os.popen(f"grep 'BAD TERMINATION' {subdir + 'no_relax.out'} | wc -l").read()) > 0
            if bad_term:
                shutil.rmtree(subdir)
                print(f"[{i}] Found bad termination at disp{j}")
                fails.append(j)
        if len(fails) > 0:
            no_fails = False
            os.chdir(pdir)
            fails = ",".join(list(map(lambda x: str(x), fails)))
            print(f"[{i}] Resubmitting failed ph jobs: {fails}")
            cmd = f'sbatch --array={fails} EXECUTABLE_BAT_DNE'
            print(f"POPEN: {cmd}"); stream = os.popen(cmd)
            print(stream.read())
            os.chdir(d)
    if no_fails:
        print("No failures in phonon force calculations.")
