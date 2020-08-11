from ___constants_names import START_BATCH_NAME, START_BATCH_PATH, BATCH_FILE_NAME, START_BATCH_OUTFILE_START, THIS_DIR

from __directory_searchers import checkPath, filesInDir
from __dirModifications import rm
from ____exit_with_error import exit_with_error

def cleanRelevantFiles():
    # Phonopy places all their files in the directory of this script. 
        # We want it in phDir, so we scane for all files that are not .py scripts or related and move them over.
        allFiles = filesInDir(THIS_DIR)
        print('Calculations for phonopy complete. Cleaning batch directory %s ...'%(START_BATCH_PATH))
        try:
            for i in allFiles:
                # The directory storing this has nothing except readme, batch file, and scripts, so all other files are phonopy
                # Move all phonopy files to the right directory
                if not (i == 'README.md' or i[-3:] == '.py' or i == START_BATCH_NAME or i == BATCH_FILE_NAME or i[:4] == START_BATCH_OUTFILE_START):
                    rm(checkPath(THIS_DIR) + i)
                    print('Removed %s from %s for cleanup.'%(i, checkPath(THIS_DIR)))
            print('Move complete.')
        except Exception as err:
            exit_with_error('Error in removing files in%s: '%(THIS_DIR) + str(err))