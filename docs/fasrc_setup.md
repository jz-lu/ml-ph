# Setting up on the FASRC Cluster

## Getting an account and learning the basics
Everything about the basics of the cluster is in the [FASRC Quickstart guide](https://docs.rc.fas.harvard.edu/kb/quickstart-guide/).

## Bashrc setup
The `.bashrc` file should be set up to include the following. The `module load` functions are required to run VASP. The path variables default the location of file access and python interpreter location. Change python version as necessary in `PYTHONPATH` (check the folder)
```
# User specific aliases and functions

module load python/3.6.3-fasrc02
module load intel/17.0.4-fasrc01
module load impi/2017.2.174-fasrc01

export PATH=/n/YOUR-HOME-DIR/YOUR-USERNAME/miniconda3/bin:$PATH

export PYTHONPATH=/n/YOUR-HOME-DIR/YOUR-USERNAME/miniconda3/lib/python3.7

# >>> conda initialize >>>
```

## Environment Setup
If running for the first time, do 
```
cd ~; wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh
$HOME/miniconda3/bin/conda create -p $HOME/ml-ph_env anaconda keras theano pygpu h5py
conda install -c conda-forge phonopy

```
to install conda (virtual environment and package manager) and build an environment with all the required packages. If any of the packages fail to install, look up the installation command on the [conda forge website](https://anaconda.org/conda-forge/repo). If prompted, let conda activate automatically on login. If `(base)` or your environment name does not show up after installation, activate conda manually with `source activate $HOME/ml-ph_env`. If you do not remember your environment name, use `conda info --envs`. You will also need another package manager, pip. Run `pip install --upgrade pip`. To get a code base started, run 
```
mkdir $HOME/codes && cd $HOME/codes && git clone https://github.com/materialsproject/pymatgen.git
cd $HOME/codes/pymatgen && pip install
```
If the last line gives an error, run `pip install pymatgen` instead. Finally, run 
```
pip install pymatgen-db pybtex bibtexparser
```
which gives all the remaining necessary packages. To deactivate conda, use `conda deactivate`. However, all work should be done within your virtual environment.

## Running on the cluster
The FASRC is split into a computing hierarchy. At the top are nodes, each of which contain a number of cores. The cluster starts in the login node, which should not be used for anything other than switching to an interactive node. Interactive nodes include `test`, as well as some Kaxiras-specific nodes (see attached file). An example command to start an interactive session is given below (we use `test` as the example node, but switch as needed). The below is quoted from the Quickstart guide.

```
srun -p test --pty --mem 500 -t 0-08:00 /bin/bash
```
"srun is like sbatch, but it runs synchronously (i.e. it does not return until the job is finished). The example starts a job on the “test” partition, with pseudo-terminal mode on (--pty), an allocation of 500 MB RAM (--mem 500), and for 6 hours (-t in D-HH:MM format). It also assumes one core on one node. The final argument is the command that you want to run. In this case you’ll just get a shell prompt on a compute host. Now you can run any normal Linux commands without taking up resources on a login node. Make sure you choose a reasonable amount of memory (--mem) for your session."

Once on the interactive session, feel free to do whatever. The run session is on a timer, but you will submit jobs using the `sbatch FILE-PATH` command, giving the batch file path in `FILE-PATH`. The batch file used to run vasp by the present tool is located in `ml-ph/ALLEGRO_ANALYZER/STATIC_BAT_DNE` and the batch file the user should call with `sbatch` is located in `ml-ph/start_job_batch/EXECUTABLE_BAT_DNE`; please do not edit either, except for the job name in the latter batch file and the number of cores/memory if more are needed. The batch file will specify the number of cores desired (make sure it is consistent with the `NUM_AVAILABLE_CORES` constant in `___constants_misc.py`) as well as the runtime. For van der Waals corrections, the time can be several hours per run, and if sampling with configuration as well, make sure to request up to 12-24 hours.