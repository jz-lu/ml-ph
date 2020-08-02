import subprocess
from ___constants_misc import unique_str

def mkdir(dirName, currPath='.'):
    if currPath[-1] != '/':
        currPath += '/'
    newDir = subprocess.run(['mkdir', currPath + dirName], capture_output=True, universal_newlines=True)
    return newDir.stderr

def move(dirName, currPath, newPath=unique_str, newName=''):
    if newPath == unique_str:
        newPath = currPath

    if currPath[-1] != '/':
        currPath += '/'
    if newPath[-1] != '/':
        newPath += '/'

    movement = subprocess.run(['mv', currPath + dirName, newPath + newName], capture_output=True, universal_newlines=True)
    return movement.stderr

def copy(dirName, currPath, newPath=unique_str, newName='', isFolder=False):
    if newPath == unique_str:
        newPath = currPath

    if currPath[-1] != '/':
        currPath += '/'
    if newPath[-1] != '/':
        newPath += '/'
    
    if isFolder:
        newCopy = subprocess.run(['cp', '-r', currPath + dirName, newPath + newName], capture_output=True, universal_newlines=True)
    else:
        newCopy = subprocess.run(['cp', currPath + dirName, newPath + newName], capture_output=True, universal_newlines=True)
    return newCopy.stderr

def rm(dirPath, isFolder=False):
    if isFolder:
        removal = subprocess.run(['rm', '-r', dirPath], capture_output=True, universal_newlines=True)
    else:
        removal = subprocess.run(['rm', dirPath], capture_output=True, universal_newlines=True)
    return removal.stderr

def cat(fileDirs, outFileName, dirOut):
    if dirOut[-1] != '/':
        dirOut += '/'

    cmd = ['cat'] + fileDirs + ['>', dirOut + outFileName]
    print(cmd)

    newFileObj = subprocess.run(cmd, capture_output=True, universal_newlines=True)
    return newFileObj.stderr

cat(['./POSCAR', './KPOINTS'], 'TESTCAT', '.')