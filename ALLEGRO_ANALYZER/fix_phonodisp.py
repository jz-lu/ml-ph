import os
from __directory_searchers import checkPath, findFilesInDir
from __dirModifications import mkdir, move
import argparse
import numpy as np

def expand_range(rg):
    try:
        idx = rg.index('-')
        start = int(rg[:idx]); end = int(rg[idx+1:])
        return np.arange(start, end+1)
    except:
        return np.array([int(rg)])

def process(dirName, supercellDim="3 3 1", Poscar_unitcell_name='POSCAR_unit'):
    dirName = checkPath(dirName)

    try:
        os.chdir(dirName)
        print(f'WD changed to "{dirName}"')
        CMD_GET_DISPLACEMENTS = ['phonopy', '-d', '--dim="{}"'.format(supercellDim), '-c', dirName + Poscar_unitcell_name] # -c flag to indicate to phonopy where unit cell POSCAR is
        CMD_GET_DISPLACEMENTS = ' '.join(CMD_GET_DISPLACEMENTS)
        print('Running "%s" to shell...'%(CMD_GET_DISPLACEMENTS))
        stream = os.popen(CMD_GET_DISPLACEMENTS)
        print(stream.read())

        poscarArray = findFilesInDir(dirName, 'POSCAR-', 'start') # Again, phonopy hard string 'POSCAR-XYZ', no need to collect
        poscarArray.sort() # Make sure displacements in order
        numPoscars = len(poscarArray)
        print('Result: {} displacement files found.'.format(numPoscars))
    except Exception as err:
        print('Phonopy preprocessing error: ' + str(err))
        return

    if numPoscars == 0:
        raise Exception("No POSCARs found")
    else:
        dispNums = []
        subdirNames = []
        # Create new directories
        print(os.popen("rm -r disp*").read())
        for i in range(numPoscars): 
            dispNums.append((poscarArray[i])[-3:]) # Gives the XYZ in POSCAR-XYZ. THIS IS A STRING!
            subdirNames.append('disp%d'%(int(dispNums[i])))
            mkdir(subdirNames[i], dirName)
            print('New subdirectory %s created.'%(checkPath(dirName + subdirNames[i])))
            move(poscarArray[i], dirName, dirName + subdirNames[i])
            assert os.path.isfile(checkPath(dirName + subdirNames[i]) + poscarArray[i])
            print(f'Moved {poscarArray[i]} to {checkPath(dirName + subdirNames[i])}')
            
        execmd = f"sbatch --array=1-{numPoscars} EXECUTABLE_BAT_DNE"
        print(f'Running "{execmd}"')
        print(os.popen(execmd).read())
        print('Job array successfully submitted')
            
        print('Total number of displacement files generated: ' + dispNums[-1])

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Resubmit erroneous-input phonopy jobs")
    parser.add_argument("sampling", type=int, help="number of shifts")
    parser.add_argument("-s", "--skip", nargs="+", help="numbers to skip", default=[])
    parser.add_argument("-d", "--dir", type=str, help="main directory", default='.')
    args = parser.parse_args()
    dir_full = checkPath(os.path.abspath(args.dir))
    
    skip = np.hstack([expand_range(x) for x in args.skip]) if len(args.skip) > 0 else []
    print(f"Skipping: {'None' if skip == [] else skip}")
    
    for i in range(args.sampling):
        if i in skip:
            continue
        process(dir_full + f'shift_{i}/analyses/phonon/')
        