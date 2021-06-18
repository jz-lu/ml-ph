# Constants for submitting jobs to a computing cluster in SLURM syntax
# must all be strings
COMPUTE_PARTITIONS = 'shared'
COMPUTE_MEM_PER_CPU = '5000'
COMPUTE_TIME = '16:00:00' # HR:MN:SC
COMPUTE_EMAIL_TO = 'jlu@college.harvard.edu'
COMPUTE_EMAIL_TYPE = 'END,FAIL' # no spaces
COMPUTE_NCPU = '2'
USE_NODE_INDICATOR = True
COMPUTE_NNODE = '1'
COMPUTE_JOBNAME = 'shifts'
PHONON_JOBNAME = 'ph'
COMPUTE_ANACONDA_ENV = 'anaconda_env' # computing requires anaconda