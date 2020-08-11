from ____exit_with_error import exit_with_error

from ___constants_names import THIS_DIR, START_BATCH_NAME, BATCH_FILE_NAME, START_BATCH_OUTFILE_START

from __directory_searchers import filesInDir
from __dirModifications import move

# For moving phonopy generated files from ths script's directory to the right one
def moveRelevantFiles(dirName):
    # Phonopy places all their files in the directory of this script. 
        # We want it in phDir, so we scane for all files that are not .py scripts or related and move them over.
        allFiles = filesInDir(THIS_DIR)
        print('Files generated by phonopy collected. Moving to the right directory ({})...'.format(dirName))
        try:
            for i in allFiles:
                # The directory storing this has nothing except readme, batch file, and scripts, so all other files are phonopy
                # Move all phonopy files to the right directory
                if not (i == 'README.md' or i[-3:] == '.py' or i == START_BATCH_NAME or i == BATCH_FILE_NAME or i[:4] == START_BATCH_OUTFILE_START):
                    move(i, THIS_DIR, dirName)
            print('Move complete.')
        except Exception as err:
            exit_with_error('Error in moving: ' + str(err))
        