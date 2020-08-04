import subprocess
import sys
from ___constants_misc import *
from ___constants_names import ROOT

# Check if a path ends in a '/' for proper use
def checkPath(dirName):
    if dirName[-1] != '/':
        dirName += '/'
    return dirName

def mkdir(dirName, currPath=ROOT):
    if currPath[-1] != '/':
        currPath += '/'
    newDir = subprocess.run(['mkdir', currPath + dirName], capture_output=True, universal_newlines=True)
    return newDir.stderr

def move(dirName, currPath, newPath=unique_str, newName=''):
    currPath = checkPath(currPath)
    if newPath == unique_str:
        newPath = currPath
    newPath = checkPath(newPath)

    movement = subprocess.run(['mv', currPath + dirName, newPath + newName], capture_output=True, universal_newlines=True)
    return movement.stderr

def copy(dirName, currPath, newPath=unique_str, newName='', isFolder=False):
    currPath = checkPath(currPath)
    newPath = checkPath(newPath)
    if newPath == unique_str:
        newPath = currPath
    
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
        sys.exit(ERR_BAD_DIR)

    return newFileObj.stdout