# Generate bash executable easily
from ___constants_compute import *
from ___constants_names import START_BATCH_NAME, ENERGIES, CMD_LINE_ARG_LIST, CODE_DIR
from __directory_searchers import checkPath
import sys, copy
from ___helpers_parsing import warn, err, is_flag, check_not_flag


args = copy.deepcopy(sys.argv)[1:]; i = 0; n = len(args)
compute_jobname = 'NoName'
compute_nnode = COMPUTE_NNODE
compute_ncpu = COMPUTE_NCPU
compute_time = COMPUTE_TIME
compute_partitions = COMPUTE_PARTITIONS
compute_mem_per_cpu = COMPUTE_MEM_PER_CPU
compute_email_type = COMPUTE_EMAIL_TYPE
compute_email_to = COMPUTE_EMAIL_TO
outdir = '.' # default to current WD
vdw = 'T'
kpts = 'GAMMA'
fname = START_BATCH_NAME
USE_NODE_INDICATOR = True
c = [ENERGIES]

while i < n:
    if not is_flag(args[i]):
        warn(f'Warning: token "{args[i]}" is out of place and will be ignored')
        i += 1; continue
    if args[i] == '-filename':
        i += 1; check_not_flag(args[i]); fname = args[i]; i +=1
    elif args[i] == '-jobname':
        i +=1; check_not_flag(args[i]); compute_jobname = args[i]; i +=1
    elif args[i] == '-n':
        i +=1; check_not_flag(args[i]); compute_ncpu = str(int(args[i])); i +=1
    elif args[i] == '-N':
        i +=1; check_not_flag(args[i]); compute_nnode = USE_NODE_INDICATOR = int(args[i]); i += 1
    elif args[i] == '-email':
        i +=1; check_not_flag(args[i]); compute_email_to = args[i]; i +=1
    elif args[i] == '-mem':
        i +=1; check_not_flag(args[i]); compute_mem_per_cpu = int(args[i]); i +=1
    elif args[i] == '-t':
        i +=1; check_not_flag(args[i]); compute_time = args[i]
        try:
            int(args[i][0:2]); int(args[i][3:5]); int(args[i][6:8])
            assert args[i][2] == args[i][5] == ':' and len(args[i]) == 8
        except:
            err('Error: specified time {args[i]} is invalid, format as HH:MM:SS')
        i += 1
    elif args[i] == '-dir':
        i +=1; check_not_flag(args[i])
        outdir = args[i]
        if args[i] in ['.', './', '..', '../'] or args[i][0] == '.':
            warn(f'Warning: specified directory "{args[i]}" may not work when running executable')
        i +=1
    elif args[i] == '-kpts':
        i +=1; check_not_flag(args[i])
        if args[i] == 'MP':
            kpts = 'MP'
        elif args[i] not in ['GAMMA', 'G', 'gamma', 'Gamma', 'Gam', 'gam', 'g']:
            warn(f'Warning: kpoints type "{args[i]}" not recognized, using Gamma-centered')
        i += 1
    elif args[i] == '-vdw':
        i +=1; check_not_flag(args[i])
        if args[i] not in ['T', 'F']:
            warn(f'Warning: expected vdW flag to be "T" or "F", got "{args[i]}", using "T"')
        i += 1
    elif args[i] == '-c':
        i += 1; check_not_flag(args[i])
        c = []
        while i < n and args[i] in CMD_LINE_ARG_LIST:
            c.append(args[i]); i += 1
        if not c or (i+1 < n and not is_flag(args[i])):
            err(f'Error: specified computations under flag "-c" invalid')
    else:
        warn(f'Warning: unknown flag "{args[i]} ignored')
        i += 1

outdir = checkPath(outdir)
c = ' '.join(c)
with open(outdir + fname, 'w') as f:
    f.write('#!/bin/bash\n')
    f.write('#SBATCH --job-name=' + compute_jobname + '\n')
    if USE_NODE_INDICATOR:
        f.write('#SBATCH -N %s\n'%(compute_nnode))
    f.write('#SBATCH -n %s\n#SBATCH -t %s\n#SBATCH -p %s\n#SBATCH --mem-per-cpu=%s\n'%(compute_ncpu, compute_time, compute_partitions, compute_mem_per_cpu))
    f.write('#SBATCH -o job_%j.out\n#SBATCH -e er_%j.err\n')
    f.write('#SBATCH --mail-type=%s\n#SBATCH --mail-user=%s\n'%(compute_email_type, compute_email_to))
    f.write('source activate $HOME/%s\n'%(COMPUTE_ANACONDA_ENV))
    f.write('WDIR="%s"\n'%(outdir))
    f.write('echo "WD: ${WDIR}"\n')
    f.write('ALLEGRO_DIR="%s"\n'%CODE_DIR)
    f.write('module list\nsource activate $HOME/%s\n'%(COMPUTE_ANACONDA_ENV))
    f.write('echo "Starting calculations..."\n')
    f.write('python3 $ALLEGRO_DIR/start.py 0 $WDIR %s %s %s\n'%(vdw, kpts, c))
    f.write('echo "Calculations complete!"\n')
print(bcolors.OKGREEN + 'Executable bash file successfully written to %s'%(outdir) + bcolors.ENDC)

