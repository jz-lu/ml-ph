from ___constants_names import START_BATCH_NAME, START_BATCH_PATH, BATCH_FILE_NAME, START_BATCH_OUTFILE_START

from __directory_searchers import checkPath, filesInDir
from __dirModifications import rm
from ____exit_with_error import exit_with_error

def cleanRelevantFiles(dirName):
    # Phonopy places all their files in the directory of this script. 
        # We want it in phDir, so we scane for all files that are not .py scripts or related and move them over.
        dirName = checkPath(dirName)
        allFiles = filesInDir(dirName)
        print('Calculations complete. Cleaning batch directory %s ...'%(START_BATCH_PATH))
        try:
            for i in allFiles:
                # The directory storing this has nothing except readme, batch file, and scripts, so all other files are phonopy
                # Move all phonopy files to the right directory
                if not (i == 'README.md' or i[-3:] == '.py' or i == START_BATCH_NAME or i == BATCH_FILE_NAME or i[:4] == START_BATCH_OUTFILE_START):
                    rm(checkPath(dirName) + i)
                    print('Removed %s from %s for cleanup.'%(i, dirName))
            print('Move complete.')
        except Exception as err:
            exit_with_error('Error in removing files in%s: '%(dirName) + str(err))