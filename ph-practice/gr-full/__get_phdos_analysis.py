from pymatgen.io.vasp.inputs import Kpoints, Poscar
import sys
import subprocess

from __directory_searchers import checkPath
from __ph_movers import moveRelevantFiles

from ___constants_names import PH_MESH_CONF_NAME
from ___constants_phonopy import SUPER_DIM, POSCAR_UNIT_NAME

# Creates the mesh.conf file necesary for phonopy dos analysis

def ph_create_mesh_conf(kpoints_mesh_obj, poscar_obj, outDir):
    outDir = checkPath(outDir)

    # Process the site symbols string from poscar
    site_symbols_str = ' '.join(poscar_obj.site_symbols)

    # Process the kpoints sampling
    kpts_mesh_str = ''
    for i in (kpoints_mesh_obj.kpts)[0]:
        kpts_mesh_str += str(i) + ' '
    kpts_mesh_str = kpts_mesh_str[:-1]
    # NOTE: this must change when not using the automatic gamma scheme to generate kpoints.

    f = open(outDir + PH_MESH_CONF_NAME, 'w')
    f.write('ATOM_NAME = ' + site_symbols_str + '\n')
    f.write('DIM = ' + SUPER_DIM + '\n')
    f.write('MP = ' + kpts_mesh_str)
    f.close()

    return

# Generate DOS raw data and plot
# Phonopy gives raw data and plot no matter what so we don't give an option here
def ph_get_dos(kpoints_mesh_obj, poscar_unitcell_obj, outDir, poscar_unit_path):
    try:
        ph_create_mesh_conf(kpoints_mesh_obj, poscar_unitcell_obj, outDir)
    except Exception as err:
        print('Error in generating mesh.conf files required for analysis fo phonon DOS.')
        sys.exit(err)

    # We need POSCAR_unit from that directory, imported from poscar_unit_path

    try:
        CMD_GET_DOS = ['phonopy', '-p', outDir + PH_MESH_CONF_NAME, '-c', poscar_unit_path, '-s']
        output_of_ph_dos_call = subprocess.run(CMD_GET_DOS, capture_output=True, universal_newlines=True)
        print(output_of_ph_dos_call)
    except Exception as err:
        sys.exit('Error in DOS analysis and output:', err)
    
    # As with before, all files generated need to be moved to this directory, by default they are in the script directory THIS_DIR
    moveRelevantFiles(outDir)

    return


