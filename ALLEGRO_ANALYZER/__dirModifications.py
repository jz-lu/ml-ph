import subprocess
from ____exit_with_error import exit_with_error
from ___constants_misc import *
from __directory_searchers import checkPath


def mkdir(dirName, currPath):
    if currPath[-1] != '/':
        currPath += '/'
    newDir = subprocess.run(['mkdir', currPath + dirName], capture_output=True, universal_newlines=True)
    if newDir.stderr != '':
        exit_with_error('Could not make new directory %s at %s. \n\n\tError: %s.\n\n'%(dirName, currPath, newDir.stderr))
    return

def move(dirName, currPath, newPath=unique_str, newName=''):
    currPath = checkPath(currPath)
    if newPath == unique_str:
        newPath = currPath
    newPath = checkPath(newPath)

    movement = subprocess.run(['mv', currPath + dirName, newPath + newName], capture_output=True, universal_newlines=True)
    if movement.stderr != '':
        exit_with_error('Could not move %s to %s. Error: %s.'%(currPath + dirName, newPath + newName, movement.stderr))
    return

def copy(dirName, currPath, newPath=unique_str, newName='', isFolder=False):
    currPath = checkPath(currPath)
    newPath = checkPath(newPath)
    if newPath == unique_str:
        newPath = currPath

    cmdarr = []
    
    if isFolder:
        cmdarr = ['cp', '-r', currPath + dirName, newPath + newName]
    else:
        cmdarr = ['cp', currPath + dirName, newPath + newName]

    newCopy = subprocess.run(cmdarr, capture_output=True, universal_newlines=True)
    print('Ran "%s"'%(' '.join(cmdarr)))
    print(newCopy.stdout)
    if newCopy.stderr != '':
        exit_with_error('Could not copy %s to %s. Error: %s.'%(currPath + dirName, newPath + newName, newCopy.stderr))
    return

def rm(dirPath, isFolder=False):
    if isFolder:
        removal = subprocess.run(['rm', '-r', dirPath], capture_output=True, universal_newlines=True)
    else:
        removal = subprocess.run(['rm', dirPath], capture_output=True, universal_newlines=True)
    if removal.stderr != '':
        exit_with_error('Could not remove %s. Error: %s.'%(dirPath, removal.stderr))
    return

def cat(fileDirs, dirOut, outFileName):
    # Takes as input an array of files (with directory paths) to concatenate, and an out directrory and new file name
    dirOut = checkPath(dirOut)

    cmd =  'cat '
    for i in fileDirs:
        cmd = cmd + i + ' '
    cmd = cmd + '> ' + dirOut + outFileName
    print(cmd)

    newFileObj = subprocess.run(cmd, shell=True, capture_output=True, universal_newlines=True)
    if newFileObj.stderr != '':
        exit_with_error(ERR_BAD_DIR)

    return newFileObj.stdout