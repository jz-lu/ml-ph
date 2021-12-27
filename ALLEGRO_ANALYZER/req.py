"""
Script to systematically search for and resubmit failed 
ml-ph VASP jobs submitted to the computing cluster.
"""
import os, shutil
import argparse

parser = argparse.ArgumentParser(description="Requeue failed jobs on cluster")
parser.add_argument("sampling", type=int, help="number of shifts")
parser.add_argument("-s", "--skip", nargs="+", help="numbers to skip", default=[])
parser.add_argument("-c", "--cfg", action="store_true", help="configuration requeue")
args = parser.parse_args()

d = os.getcwd()
skip = list(map(lambda x: int(x), args.skip))
if len(skip) > 0:
    print(f"Skipping: {skip}")
N = args.sampling
if args.cfg:
    fails = []
    for i in range(N):
        if i in skip:
            continue
        s1 = os.popen(f"grep 'not available' shift_{i}/relaxation/relax.err | wc -l").read()
        s2 = os.popen(f"grep 'rror' shift_{i}/relaxation/relax.err | wc -l").read()
        if int(s1) > 0 or int(s2) > 0:
            if os.path.isfile(f'shift_{i}/relaxation/CONTCAR'):
                shutil.copy(f'shift_{i}/relaxation/CONTCAR', f'shift_{i}/POSCAR')
            shutil.rmtree(f'shift_{i}/relaxation/')
            shutil.rmtree(f'shift_{i}/analyses/')
            fails.append(i)
    if len(fails) > 0:
        fails = ",".join(list(map(lambda x: str(x), fails)))
        print(f"Requeueing fails: {fails}")
        stream = os.popen(f"sbatch --array={fails} EXECUTABLE_BAT_DNE")
        print(stream.read())
    else:
        print("No failures detected in configuration relaxations.")
    
else:
    no_fails = True
    for i in range(N):
        if i % 10 == 0:
            print(f"Analysis finished up to shift {i}")
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
            other_err = int(os.popen(f"grep 'rror' {subdir + 'no_relax.err'} | wc -l").read()) > 0
            unknown_failure = int(os.popen(f"grep 'fail' {subdir + 'no_relax.err'} | wc -l").read()) > 0
            never_finished = int(os.popen(f"tail {subdir + 'no_relax.out'} -n 1 | grep DAV | wc -l").read()) > 0
            if bad_term or other_err or unknown_failure or never_finished:
                shutil.rmtree(subdir)
                if bad_term:
                    print(f"[{i}] Found bad termination at disp{j}")
                elif other_err or unknown_failure:
                    print(f"[{i}] Found unknown error at disp{j}")
                elif never_finished:
                    print(f"[{i}] Found job never finished at disp{j}")
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
        print("No failures detected in phonon force calculations.")
