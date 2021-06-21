import os, fnmatch
from ____exit_with_error import exit_with_error
from ___constants_misc import ERR_INVALID_FINDDIR_PARAM

# Not currently used for anything.
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

# Find all files in a given directory (no subdirectory search, use find() for that)
def filesInDir(dirName):
    dirName = checkPath(dirName)
    files = [f for f in os.listdir(dirName) if os.path.isfile(os.path.join(dirName, f))]
    return files

# Find all files in a directory. Specify type as 'start', 'end', or 'exact' to get search that starts with, ends with, or is exactly, fileName.
#   Return files in alphanumerically ordered list.
# def findFilesInDir(dirName, fileName, searchType='exact'):
#     dirName = checkPath(dirName)
#     print("Searching for matches '%s' of file '%s' in path %s with files:"%(searchType, fileName, dirName), 
#         os.listdir(dirName))
#     arr = []
    
#     if searchType == 'exact':
#         for f in os.listdir(dirName):
#             if f == fileName:
#                 arr.append(f)

#     elif searchType == 'start':
#         for f in os.listdir(dirName):
#             if f.startswith(fileName):
#                 arr.append(f)
    
#     elif searchType == 'end':
#         for f in os.listdir(dirName):
#             if f.endswith(fileName):
#                 arr.append(f)
    
#     else:
#         exit_with_error(ERR_INVALID_FINDDIR_PARAM)

#     arr.sort()
#     print("Search done, returning list:", arr)
#     return arr

def findFilesInDir(dirName, fileName, searchType='exact'):
    dirName = checkPath(dirName)
    print("Searching for matches '%s' of filename'%s' in directory %s with list:"%(searchType, fileName, dirName), 
        os.listdir(dirName))
    arr = []
    
    if searchType == 'exact':
        arr += list(filter(lambda x: x == fileName and os.path.isfile(x), os.listdir()))
    elif searchType == 'start':
        arr += list(filter(lambda x: x.startswith(fileName) and os.path.isfile(x), os.listdir()))
    elif searchType == 'end':
        arr += list(filter(lambda x: x.endswith(fileName) and os.path.isfile(x), os.listdir()))
    else:
        exit_with_error(ERR_INVALID_FINDDIR_PARAM)

    arr.sort()
    print("Search done, returning list:", arr)
    return arr

def findDirsinDir(dirName, subName, searchType='exact'):
    dirName = checkPath(dirName)
    print("Searching for matches '%s' of subdirectory '%s' in directory %s with list:"%(searchType, subName, dirName), 
        os.listdir(dirName))
    arr = []
    
    if searchType == 'exact':
        arr += list(filter(lambda x: x == subName and os.path.isdir(x), os.listdir()))
    elif searchType == 'start':
        arr += list(filter(lambda x: x.startswith(subName) and os.path.isdir(x), os.listdir()))
    elif searchType == 'end':
        arr += list(filter(lambda x: x.endswith(subName) and os.path.isdir(x), os.listdir()))
    else:
        exit_with_error(ERR_INVALID_FINDDIR_PARAM)

    arr.sort()
    print("Search done, returning list:", arr)
    return arr
