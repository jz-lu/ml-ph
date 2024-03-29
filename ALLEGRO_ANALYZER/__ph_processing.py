from ____exit_with_error import exit_with_error
from __directory_searchers import find # Will let us look for the different POSCAR-XYZ displacement files
from __dirModifications import move, copy, rm, mkdir
from __directory_searchers import checkPath, findFilesInDir, filesInDir
from __run_vasp import run_vasp
from __input_modifiers import modifyIncar
from __ph_movers import moveRelevantFiles
from ___constants_phonopy import *
from ___constants_names import *
from ___constants_misc import ERR_PH_FORCE_SETS, ERR_PH_FORCE_SETS_NOT_FOUND
from ___constants_vasp import *
from ___constants_compute import *
from make_execscr import build_bash_exe
import subprocess
import os
from pymatgen.io.vasp.inputs import Poscar, VaspInput, Potcar # pylint: disable=import-error

# Builds and enqueues batch file for computing phonon displacements. BASE_ROOT should specify the path to the given displacment
def compute_displacements(ROOT, user_input_settings, ndisp):
    assert os.path.isdir(ROOT), f"Directory {ROOT} does not exist"
    ROOT = checkPath(ROOT)
    print(f'Writing a job-array bash executable to {ROOT}')
    vdw = 'T' if user_input_settings.do_vdW in ['T', True] else 'F'
    subname = user_input_settings.passname()
    large_sc = (user_input_settings.get_super_dim() > 3) # larger supercells need more CPUs
    monster = (large_sc and (vdw == 'T')) # hardest type of calculation to converge
    if user_input_settings.get_type_flag() == TYPE_RELAX_CONFIG:
        if user_input_settings.sampling_is_diagonal():
            subname += '-d'
        else:
            subname += '-c'
    if vdw == 'T':
        subname += '-v'
    kpts = 'GAMMA' if user_input_settings.kpoints_is_gamma_centered else 'MP'
    exepath = build_bash_exe(calc_type='norelax', outdir=ROOT, wdir=ROOT+PHDISP_STATIC_NAME, 
                   calc_list=[ENERGIES], compute_jobname=PHONON_JOBNAME+subname, 
                   compute_mem_per_cpu='8000', 
                   vdw=vdw, kpts=kpts, 
                   compute_time='24:00:00' if monster else COMPUTE_TIME, 
                   as_arr=True, compute_ncpu='72' if large_sc else '48',
                   USE_NODE_INDICATOR=False if large_sc else True)
    print('Executable built.')
    runcmd = 'sbatch --array=1-%d'%(ndisp) + ' ' + exepath
    print('Running %s...'%runcmd)
    stream = os.popen(runcmd)
    print(stream.read())
    return None

# Compute all of the forces for a given displacement-generated phononic supercell via VASP.
def compute_vasp_ph_forces(index, dispNum, dirName, subdirName, disp_poscar_name, vaspObj):
    try:
        # Perform the VASP calculation on this 
        # NOTE: no chgcar, the charge densities would be inaccurate now that we have displacements
        displaced_poscar = Poscar.from_file(dirName + disp_poscar_name)
        ph_vasp_obj = VaspInput(vaspObj['INCAR'], vaspObj['KPOINTS'], displaced_poscar, vaspObj['POTCAR'])
        print('[disp %d] Starting VASP nonrelaxation forces calculation for displacement '%(index+1) + dispNum + '...')
        print('Calculations to be stored in %s'%(dirName + subdirName))
        run_vasp(ph_vasp_obj, dirName + subdirName, run_type='phonon')
        print('VASP calculation complete.')
    except Exception as err:
        print('An error occurred while processing displacement POSCAR file {}'.format(dirName + disp_poscar_name))
        exit_with_error('Error in preprocessing phonopy (parsing displacement files and running VASP force calculations): ' + str(err))

# Function to handle all preprocessing of phonopy before the force calulation.
def ph_preprocess(dirName, vaspObj, supercellDim=SUPER_DIM_STR, Poscar_unitcell_name=POSCAR_UNIT_NAME, 
                  onejob=True, user_input_settings=None):
    print(PHONOPY_DISP_MSG)
    dirName = checkPath(dirName)
    assert (user_input_settings is not None) ^ onejob

    try:
        os.chdir(dirName)
        print(f'WD changed to "{dirName}"')
        CMD_GET_DISPLACEMENTS = ['phonopy', '-d', '--dim="{}"'.format(supercellDim), '-c', dirName + Poscar_unitcell_name] # -c flag to indicate to phonopy where unit cell POSCAR is
        CMD_GET_DISPLACEMENTS = ' '.join(CMD_GET_DISPLACEMENTS)
        print('Running "%s" to shell...'%(CMD_GET_DISPLACEMENTS))
        stream = os.popen(CMD_GET_DISPLACEMENTS)
        print(stream.read())

        # Phonopy places all their files in the directory of this script. 
        # We want it in phDir, so we scane for all files that are not .py scripts or related and move them over.
        # moveRelevantFiles(dirName)

        # ...Except you need phonopy_disp.yaml in the current directory for FORCE_SETS. We'll clean it after it's generated.
        # copy(PH_DISP_YAML_NAME, dirName, THIS_DIR)
        # print('From ph_processing module: Copy of %s complete.'%(PH_DISP_YAML_NAME))
        # print('Verification of existence in %s :'%(THIS_DIR), findFilesInDir(THIS_DIR, PH_DISP_YAML_NAME))

        # Check for the number of POSCAR-XYZ's made and organize them accordingly.
        poscarArray = findFilesInDir(dirName, POSCAR_PH_NAMEPRE, 'start') # Again, phonopy hard string 'POSCAR-XYZ', no need to collect
        poscarArray.sort() # Make sure displacements in order
        numPoscars = len(poscarArray)
        print('Result: {} displacement files found.'.format(numPoscars))
    except Exception as err:
        exit_with_error('Phonopy preprocessing error: ' + str(err))

    if numPoscars == 0:
        exit_with_error(PHONOPY_DISP_ERR_1)
    elif numPoscars > 999:
        exit_with_error(PHONOPY_DISP_ERR_2)
    else:
        dispNums = []
        subdirNames = []
        # Create new directories
        for i in range(numPoscars): 
            dispNums.append((poscarArray[i])[-3:]) # Gives the XYZ in POSCAR-XYZ. THIS IS A STRING!
            subdirNames.append(PHDISP_DIR_NAME%(int(dispNums[i])))
            mkdir(subdirNames[i], dirName)
            print('New subdirectory %s created.'%(checkPath(dirName + subdirNames[i])))
            if not onejob:
                move(poscarArray[i], dirName, dirName + subdirNames[i])
                assert os.path.isfile(checkPath(dirName + subdirNames[i]) + poscarArray[i])
            print(f'Moved {poscarArray[i]} to {checkPath(dirName + subdirNames[i])}')

        if onejob:
            print('Computing forces via VASP in this job')
            print("Submitting the following INCAR to VASP for phonon calculation:\n", vaspObj['INCAR'])
            for i in range(numPoscars):
                compute_vasp_ph_forces(i, dispNums[i], dirName, subdirNames[i], poscarArray[i], vaspObj)
        else:
            print('Submitting a new job array to compute forces for each displacement...')
            compute_displacements(dirName, user_input_settings, numPoscars)
            print('Job array successfully submitted')
            
        print('Total number of displacement files generated: ' + dispNums[-1])

    return dispNums[-1] # Returns last displacement number so we can get force sets.

# Takes as input dispNum, the output of phPreProcess
def ph_generate_forcesets(dirName, dispNum, path_pad=None, dump=False):
    assert int(dispNum) >= 1, f"Invalid number of displacements {int(dispNum)}: must be at least 1"
    dispNum = str(dispNum); dirName = checkPath(dirName); os.chdir(dirName) # Move to the new directory
    if dump:
        print('CHANGED WD to %s for storing phonopy calculations.'%(dirName))
        print('Files/directories in CWD:', os.listdir())
    try:
        pad = checkPath(path_pad) if path_pad else ''
        cmd = 'phonopy -f %s{1..%s}/%s%s -c %s'%(dirName + PHDISP_STATIC_NAME, dispNum, pad, VASP_RUN_XML_NAME, POSCAR_UNIT_NAME)
        if dump:
            print('Generating force sets now...')
            print('Running command to shell: ' + cmd)
        phForce_output = subprocess.run(cmd, shell=True, capture_output=True, universal_newlines=True)
        print(phForce_output.stdout)
        print(phForce_output.stderr)
    except Exception as err:
        print(ERR_PH_FORCE_SETS)
        print('Command log: {}'.format(phForce_output.stderr))
        exit_with_error('Error: ' + str(err))

    # If force sets created then proceed, if not exit with error. Use find function to find out.
    if not os.path.isfile(dirName + PH_FORCE_SETS_NAME):
        exit_with_error(ERR_PH_FORCE_SETS_NOT_FOUND)
    print("FORCE_SETS put in directory " + dirName)
    return 
    
# Get a neat function that uses the above functions to handle all processing
def ph_prepare_for_analysis(rootDirName, incar_selfcon, kpoints_mesh_nonrelax, poscar_relaxed, potcar):
    print('Phonopy analysis starting...')
    rootDirName = checkPath(rootDirName)
    mkdir(PHONOPY_DIR_NAME, rootDirName)
    DIR_PHONOPY = checkPath(rootDirName + PHONOPY_DIR_NAME)
    os.chdir(DIR_PHONOPY) # Move to the new directory
    print('os.chdir() CALLED to cd to %s for storing phonopy calculations.'%(DIR_PHONOPY))

    # Write out the poscar into the phonopy directory we made, so we can use it for calling phonopy
    try:
        poscar_relaxed.write_file(DIR_PHONOPY + POSCAR_UNIT_NAME)
    except Exception as err:
        print('Error in writing the %s file to %s.'%(POSCAR_UNIT_NAME, DIR_PHONOPY))
        exit_with_error(err)

    # Due to the displacement invalidating charge densities, it is important that we use default charge densities to start
    # i.e. ICHARG = default = 2
    incar_selfcon_initChg_noWriteChg = modifyIncar(incar_selfcon, addArr=[('ICHARG', ICHARG['default']), ('LCHARG', LCHARG['no_write_charge'])])

    # Note that in this object the poscar could be any valid poscar; we'll replace it in preprocessing by displacement poscar
    ph_preprocess_vasp_obj = VaspInput(incar_selfcon_initChg_noWriteChg, kpoints_mesh_nonrelax, poscar_relaxed, potcar)

    # Run preprocessing
    print('Sending vasp object for phonopy-specific force calculations to phonopy preprocessing module...')
    #print('[DEBUGMSG] Phonopy vasp object:', ph_preprocess_vasp_obj)
    lastDispNum = ph_preprocess(DIR_PHONOPY, ph_preprocess_vasp_obj) # returns a string with largest XYZ in POSCAR-XYZ for use in force sets generator
    
    # Generate force sets file
    ph_generate_forcesets(DIR_PHONOPY, lastDispNum)

    return DIR_PHONOPY

    