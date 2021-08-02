# Automated bash executable generator for Harvard FASRC computing cluster
from ___constants_compute import (
    COMPUTE_PARTITIONS, 
    COMPUTE_MEM_PER_CPU, 
    COMPUTE_TIME, 
    COMPUTE_EMAIL_TO, 
    COMPUTE_EMAIL_TYPE, 
    COMPUTE_NCPU, 
    USE_NODE_INDICATOR, 
    COMPUTE_NNODE, 
    COMPUTE_JOBNAME, 
    COMPUTE_ANACONDA_ENV
)
from ___constants_names import (
    START_BATCH_NAME, ENERGIES, CMD_LINE_ARG_LIST, CODE_DIR, 
    TYPE_STRS, CFG_STRS, 
    TYPE_RELAX_BASIC, TYPE_RELAX_CONFIG, TYPE_TWISTED_CONFIG, TYPE_NORELAX_BASIC
)
from __directory_searchers import checkPath
import sys, copy
from ___helpers_parsing import greet, succ, warn, err, is_flag, check_not_flag


def build_bash_exe(calc_type='basic', outdir='.', wdir=None, calc_list=[ENERGIES], 
                   compute_jobname='NoName', compute_nnode=COMPUTE_NNODE, compute_ncpu=COMPUTE_NCPU, 
                   compute_time=COMPUTE_TIME, compute_partitions=COMPUTE_PARTITIONS, 
                   compute_mem_per_cpu=COMPUTE_MEM_PER_CPU, compute_email_type=COMPUTE_EMAIL_TYPE, 
                   compute_email_to=COMPUTE_EMAIL_TO, vdw=False, kpts='GAMMA', fname=START_BATCH_NAME,
                   USE_NODE_INDICATOR=True, as_arr=False, twist=None, sampling='low', passname='', 
                   pass_idx=False):
    assert calc_type in TYPE_STRS, f"Unknown calculation type {calc_type}"
    assert isinstance(calc_list, list), "Calculation list must be of type list"
    assert vdw in ['T', 'F', True, False], "vdw parameter must be either 'T' or 'F'"
    assert kpts in ['GAMMA', 'G', 'Gamma', 'Gam', 'g', 'M', 'MP', True, False], "k-points parameter must be either 'GAMMA' or 'MP'"
    assert sampling in ['low', 'high'], "sampling must be low or high"
    calc_list = ' '.join(calc_list)
    greet(f"Building new bash exe file of type {calc_type} (calc list={calc_list})...")
    outdir = checkPath(outdir); exepath = outdir + fname
    wdir = outdir if wdir is None else wdir
    kpts = '' if kpts in ['GAMMA', 'G', 'Gamma', 'Gam', 'g', True] else '--mp'
    vdw = '--vdw' if vdw in ['T', True, 'True'] else ''
    twist = '' if twist is None else f'--twist {twist}'
    sampling = f'-s {sampling}' if calc_type in CFG_STRS else ''
    if pass_idx:
        passname += ('' if passname == '' else '-') + '${SLURM_ARRAY_TASK_ID}'
    passname = '-n ' + passname

    with open(exepath, 'w') as f:
        f.write('#!/bin/bash\n')
        if as_arr:
            f.write('#SBATCH -J %s\n'%(compute_jobname))
        else:
            f.write('#SBATCH --job-name=' + compute_jobname + '\n')
        if USE_NODE_INDICATOR:
            f.write('#SBATCH -N %s\n'%(compute_nnode))
        f.write('#SBATCH -n %s\n#SBATCH -t %s\n#SBATCH -p %s\n#SBATCH --mem-per-cpu=%s\n'%(compute_ncpu, compute_time, compute_partitions, compute_mem_per_cpu))
        if as_arr:
            f.write('#SBATCH -o %s_%%a_%%A.out\n#SBATCH -e er_%s_%%a_%%A.err\n'%(compute_jobname, compute_jobname))
        else:
            f.write('#SBATCH -o %s_%%j.out\n#SBATCH -e er_%s_%%j.err\n'%(compute_jobname, compute_jobname))
        f.write('#SBATCH --mail-type=%s\n#SBATCH --mail-user=%s\n'%(compute_email_type, compute_email_to))
        f.write('source activate $HOME/%s\n'%(COMPUTE_ANACONDA_ENV))
        if as_arr:
            f.write('WDIR="%s${SLURM_ARRAY_TASK_ID}"\n'%(wdir))
        else:
            f.write('WDIR="%s"\n'%(outdir))
        f.write('echo "WD: ${WDIR}"\n')
        f.write('ALLEGRO_DIR="%s"\n'%CODE_DIR)
        f.write('module load julia\nmodule list\nsource activate $HOME/%s\n'%(COMPUTE_ANACONDA_ENV))
        f.write('echo "Starting calculations..."\n')
        f.write(f'python3 $ALLEGRO_DIR/start.py -t {calc_type} {twist} {sampling} -d $WDIR {passname} {vdw} {kpts} {calc_list}\n')
        f.write('echo "Calculations complete!"\n')
    succ('Executable bash file successfully written to %s'%(exepath))
    return exepath

if __name__ == '__main__':
    args = copy.deepcopy(sys.argv)[1:]; i = 0; n = len(args)
    compute_jobname = COMPUTE_JOBNAME
    compute_nnode = COMPUTE_NNODE
    compute_ncpu = COMPUTE_NCPU
    compute_time = COMPUTE_TIME
    compute_partitions = COMPUTE_PARTITIONS
    compute_mem_per_cpu = COMPUTE_MEM_PER_CPU
    compute_email_type = COMPUTE_EMAIL_TYPE
    compute_email_to = COMPUTE_EMAIL_TO
    outdir = '.' # default to current WD
    vdw = ''
    sampling = 'low'
    kpts = ''
    fname = START_BATCH_NAME
    USE_NODE_INDICATOR = True
    calc_list = [ENERGIES]
    calc_type = 'basic'
    twist = None

    while i < n:
        if not is_flag(args[i]):
            warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
            i += 1; continue
        if args[i] == '-filename':
            i += 1; check_not_flag(args[i]); fname = args[i]; i += 1
        elif args[i] == '-type':
            i += 1; check_not_flag(args[i]); assert args[i] in TYPE_STRS; calc_type = args[i]; i += 1
        elif args[i] == '-jobname':
            i += 1; check_not_flag(args[i]); compute_jobname = args[i]; i += 1
        elif args[i] == '-n':
            i += 1; check_not_flag(args[i]); compute_ncpu = str(int(args[i])); i += 1
        elif args[i] == '-N':
            i += 1; check_not_flag(args[i]); compute_nnode = USE_NODE_INDICATOR = int(args[i]); i += 1
        elif args[i] == '-email':
            i += 1; check_not_flag(args[i]); compute_email_to = args[i]; i += 1
        elif args[i] == '-mem':
            i += 1; check_not_flag(args[i]); compute_mem_per_cpu = int(args[i]); i += 1
        elif args[i] == '-t':
            i += 1; check_not_flag(args[i]); compute_time = args[i]
            try:
                int(args[i][0:2]); int(args[i][3:5]); int(args[i][6:8])
                assert args[i][2] == args[i][5] == ':' and len(args[i]) == 8
            except:
                err('Error: specified time {args[i]} is invalid, format as HH:MM:SS')
            i += 1
        elif args[i] == '-dir':
            i += 1; check_not_flag(args[i])
            outdir = args[i]
            if args[i] in ['.', './', '..', '../'] or args[i][0] == '.':
                warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
            i += 1
        elif args[i] == '-s':
            i += 1; check_not_flag(args[i]); assert args[i] in ['low', 'high']
            sampling = args[i]; i += 1
        elif args[i] == '-kpts':
            i += 1; check_not_flag(args[i])
            if args[i] == 'MP':
                kpts = '--mp'
            elif args[i] not in ['GAMMA', 'G', 'gamma', 'Gamma', 'Gam', 'gam', 'g']:
                warn(f'Warning: kpoints type "{args[i]}" not recognized, using Gamma-centered')
            i += 1
        elif args[i] == '-vdw':
            i += 1; check_not_flag(args[i])
            vdw = '--vdw' if args[i] == 'T' else ''
            if args[i] not in ['T', 'F']:
                warn(f'Warning: expected vdW flag to be "T" or "F", got "{args[i]}", using "F"')
            i += 1
        elif args[i] == '-tw':
            i += 1; check_not_flag(args[i]); assert 0 < args[i] < 180
            twist = float(args[i]); i += 1
        elif args[i] == '-c':
            i += 1; check_not_flag(args[i])
            calc_list = []
            while i < n and args[i] in CMD_LINE_ARG_LIST:
                calc_list.append(args[i]); i += 1
            if not calc_list or (i+1 < n and not is_flag(args[i])):
                err(f'Error: specified computations under flag "-c" invalid')
        else:
            warn(f'Warning: unknown flag "{args[i]} ignored')
            i += 1

    build_bash_exe(calc_type=calc_type, outdir=outdir, calc_list=calc_list, compute_jobname=compute_jobname,
                    compute_nnode=compute_nnode, compute_ncpu=compute_ncpu, 
                    compute_time=compute_time, compute_partitions=compute_partitions, 
                    compute_mem_per_cpu=compute_mem_per_cpu, compute_email_type=compute_email_type, 
                    compute_email_to=compute_email_to, kpts=kpts, vdw=vdw, fname=fname,
                    USE_NODE_INDICATOR=USE_NODE_INDICATOR, twist=twist, sampling=sampling)

