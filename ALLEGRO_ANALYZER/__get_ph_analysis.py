from pymatgen.io.vasp.inputs import Kpoints, Poscar
import subprocess
import numpy as np

from ____exit_with_error import exit_with_error

from ___constants_names import PH_MESH_CONF_NAME, PH_BAND_CONF_NAME, KPOINTS_LINE_NAME
from ___constants_phonopy import SUPER_DIM, POSCAR_UNIT_NAME
from ___constants_misc import ERR_PH_CANNOT_GEN_MESHCONF
from ___constants_phonopy import PHONOPY_NUM_LINE_INTS

from __directory_searchers import checkPath
from __ph_movers import moveRelevantFiles

# Creates the mesh.conf file necessary for phonopy dos analysis
def ph_create_mesh_conf(kpoints_mesh_obj, poscar_obj, outDir):
    outDir = checkPath(outDir)

    # Process the site symbols string from poscar
    print('Building mesh.conf file...')
    site_symbols_str = ' '.join(poscar_obj.site_symbols)

    # Process the kpoints sampling
    kpts_mesh_str = ''
    for i in (kpoints_mesh_obj.kpts)[0]:
        kpts_mesh_str += str(i) + ' '
    kpts_mesh_str = kpts_mesh_str[:-1]

    if kpts_mesh_str == '0 0 0' or kpts_mesh_str == '0.0 0.0 0.0':
        exit_with_error(ERR_PH_CANNOT_GEN_MESHCONF)
    # NOTE: this must change when not using the automatic gamma scheme to generate kpoints.

    f = open(outDir + PH_MESH_CONF_NAME, 'w')
    f.write('ATOM_NAME = ' + site_symbols_str + '\n')
    f.write('DIM = ' + SUPER_DIM + '\n')
    f.write('MESH = ' + kpts_mesh_str + '\n')

    # If the kpoints is Gamma-centered, then we set gamma centered to TRUE
    if str(kpoints_mesh_obj.style).lower()[0] == 'g':
        f.write('GAMMA_CENTER = .TRUE.' + '\n')

    # Any shift we need can also be added here now
    kshift = ''
    for i in kpoints_mesh_obj.kpts_shift:
        kshift += str(i) + ' '
    kshift = kshift[:-1]
    if kshift != '0 0 0': # i.e. it's a nonzero shift
        f.write('MP_SHIFT = ' + kshift)

    f.close()

    print('mesh.conf built. Written to %s'%(outDir + PH_MESH_CONF_NAME))

    return

# The next two functions are necessary to create band.conf
# Get the path string for phonopy from LINE_KPOINTS
def ph_get_band_path(kpoints_line_obj):
    print('Getting phononpy band path from Kpoints line object for band.conf...')
    # Fancy way of removing duplicate rows, which VASP needs but phonopy does not want
    # NOTE: this function doesn't work well for complicated paths, but 2D materials are probably fine. See github's homepage README for more info.
    new_array = kpoints_line_obj.kpts
    np_path = np.unique(new_array, axis=0)
    path = []
    for i in np_path:
        path.append(list(i))
    path.append(path[0]) # We do want the last row to be the same as the first row, to get a closed path!

    # To join it all together into a path string 
    # formatted by single space between path coordinates in a point and oduble space between each path point
    for i in range(0, len(path)):
        for j in range(0, len(path[i])):
            path[i][j] = str(path[i][j])

    for i in range(0, len(path)):
        path[i] = ' '.join(path[i])
    
    path = '  '.join(path)
    print('Path made: ' + path)

    return path

# Get the labels for the path
def ph_get_band_path_labels(kpoints_line_obj):
    print('Getting phonopy band path labels for band.conf...')
    # We really just need to get rid of any labels that don't 
    # NOTE: this function doesn't work well for complicated paths, but 2D materials are probably fine. See github's homepage README for more info.
    labelArr = list(dict.fromkeys(kpoints_line_obj.labels))
    gammaIndex = None

    # phonopy expects a special input style for Gamma different from vasp
    for i in range(0, len(labelArr)):
        if labelArr[i] == '\\Gamma':
            gammaIndex = i
    labelArr[gammaIndex] = '$\Gamma$'

    # like with get_band_path, we append to the end the first label to close the path
    labelArr.append(labelArr[0])

    print('Path labels created: ' + ' '.join(labelArr))

    # Mash it together into a string
    return ' '.join(labelArr)

# Creates the band.cond file necessary for phonopy bana analysis
def ph_create_band_conf(kpoints_line_obj, poscar_obj, outDir):
    outDir = checkPath(outDir)

    # Process the site symbols string from poscar
    site_symbols_str = ' '.join(poscar_obj.site_symbols)

    # Write it out
    print('Building band.conf file...')
    f = open(outDir + PH_BAND_CONF_NAME, 'w')
    f.write('ATOM_NAME = ' + site_symbols_str + '\n')
    f.write('DIM = ' + SUPER_DIM + '\n')
    f.write('BAND = ' + ph_get_band_path(kpoints_line_obj) + '\n')
    # TODO: uncomment this when you figure it out
    # f.write('BAND_LABELS = ' + ph_get_band_path_labels(kpoints_line_obj) + '\n')
    if PHONOPY_NUM_LINE_INTS > 51:
        f.write('BAND_POINTS = ' + str(PHONOPY_NUM_LINE_INTS) + '\n')
    f.close()

    print('band.conf built. Written to %s'%(outDir + PH_BAND_CONF_NAME))

    return

# Generate DOS raw data and plot
# Phonopy gives raw data and plot no matter what so we don't give an option here
# NOTE: the unit poscar object is to get the name for the .conf file and the path is to run the subprocess in the actual generation
def ph_get_dos(kpoints_mesh_obj, poscar_unitcell_obj, outDir, poscar_unit_path):
    outDir = checkPath(outDir)
    try:
        ph_create_mesh_conf(kpoints_mesh_obj, poscar_unitcell_obj, outDir)
    except Exception as err:
        print('Error in generating mesh.conf files required for analysis fo phonon DOS.')
        exit_with_error(err)

    # We need POSCAR_unit from that directory, imported from poscar_unit_path

    try:
        CMD_GET_DOS = ['phonopy', '-p', outDir + PH_MESH_CONF_NAME, '-c', poscar_unit_path, '-s']
        print('Running command to shell:', ' '.join(CMD_GET_DOS))
        output_of_ph_dos_call = subprocess.run(CMD_GET_DOS, capture_output=True, universal_newlines=True)
        print(output_of_ph_dos_call.stdout)
    except Exception as err:
        print(output_of_ph_dos_call.stderr)
        exit_with_error('Error in DOS analysis and output: ' + str(err))
    
    # As with before, all files generated need to be moved to this directory, by default they are in the script directory THIS_DIR
    moveRelevantFiles(outDir)

    return

# Get band structure raw data and plot, as with DOS it's a nonnegotiable package deal
def ph_get_band(kpoints_line_obj, poscar_unitcell_obj, outDir, poscar_unit_path):
    outDir = checkPath(outDir)
    try:
        ph_create_band_conf(kpoints_line_obj, poscar_unitcell_obj, outDir)
    except Exception as err:
        print('Error in generating mesh.conf files required for analysis fo phonon DOS.')
        exit_with_error(err)

    # We need POSCAR_unit from that directory, imported from poscar_unit_path

    try:
        CMD_GET_BAND = ['phonopy', '-p', outDir + PH_BAND_CONF_NAME, '-c', poscar_unit_path, '-s']
        print('Running command to shell:', ' '.join(CMD_GET_BAND))
        output_of_ph_band_call = subprocess.run(CMD_GET_BAND, capture_output=True, universal_newlines=True)
        print(output_of_ph_band_call.stdout)
    except Exception as err:
        print(output_of_ph_band_call.stderr)
        exit_with_error('Error in phonon band analysis and output: ' + str(err))
    
    # As with before, all files generated need to be moved to this directory, by default they are in the script directory THIS_DIR
    moveRelevantFiles(outDir)

    return

# TEST CASES
# k = Kpoints.from_file('LINE_KPOINTS')
# kmesh = Kpoints.from_file('KPOINTS')
# p = Poscar.from_file('POSCAR_unit')
# # ph_get_dos(kmesh, p, '.', './POSCAR_unit')
# ph_get_band(k, p, '.', './POSCAR_unit')