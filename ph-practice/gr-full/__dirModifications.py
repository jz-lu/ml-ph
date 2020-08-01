import subprocess
_unique_str = 'aUniqueStringThatIsVeryUnlikelyToBeReproducedRandomly_q3452435234'

def mkdir(dirName, currPath='.'):
    if currPath[-1] != '/':
        currPath += '/'
    newDir = subprocess.run(['mkdir', currPath + dirName], capture_output=True, universal_newlines=True)
    return newDir.stderr

def move(dirName, currPath, newPath=_unique_str, newName=''):
    if newPath == _unique_str:
        newPath = currPath

    if currPath[-1] != '/':
        currPath += '/'
    if newPath[-1] != '/':
        newPath += '/'

    movement = subprocess.run(['mv', currPath + dirName, newPath + newName], capture_output=True, universal_newlines=True)
    return movement.stderr

def copy(dirName, currPath, newPath=_unique_str, newName='', isFolder=False):
    if newPath == _unique_str:
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

def cat(name1, name2, outFileName, dir1, dir2=_unique_str, dirOut=_unique_str):
    if dir1[-1] != '/':
        dir1 += '/'
    if dir2 == _unique_str:
        dir2 = dir1
    if dirOut=_unique_str:
        dirOut = dir1
    newFileObj = subprocess.run(['cat', dir1 + name1, dir2 + name2, '>', dirOut + outFileName], capture_output=True, universal_newlines=True)
    return newFileObj.stderr