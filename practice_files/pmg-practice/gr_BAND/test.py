# A find function we'll need later for finding how many POSCARs there are
import os, fnmatch
def find(pattern, path):
    result = []
    for root, dirs, files in os.walk(path):
        for name in files:
            if fnmatch.fnmatch(name, pattern):
                result.append(name)
    return result

RELATIVE_PH_DIR = '/Users/jonathanlu/Documents/ml-ph/pmg-practice/gr_BAND'
poscarArray = find('POSCAR-*', RELATIVE_PH_DIR)
print(poscarArray)