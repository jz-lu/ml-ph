# A function we'll need later for finding how many POSCARs there are
import os, fnmatch
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(name)
    return result

# Check if a path ends in a '/' for proper use
def checkPath(dirName):
    if dirName[-1] != '/':
        dirName += '/'
    return dirName