from pymatgen.io.vasp.outputs import Outcar
from ____exit_with_error import exit_with_error
from ___constants_names import TOT_ENERGIES_NAME, OUTCAR_NAME, TOTAL_ENER_DIR_NAME
from __directory_searchers import checkPath
from __dirModifications import mkdir
from ___constants_misc import ERR_ENER_WRITE_FAIL

# Input calculation output object, output tuple (tot_energy, fermi_energy)
def get_energy_analysis(dirName, outcar, writeOut=False):
    dirName = checkPath(dirName)
    energyPair = (outcar.final_energy, outcar.efermi)
    try: 
        energyFile = open(dirName + TOT_ENERGIES_NAME, 'w')

        str1 = 'Total energy (eV): ' + str(outcar.final_energy)
        str2 = 'Fermi energy (eV): ' + str(outcar.efermi)

        strings = [str1, '\n', str2]

        energyFile.writelines(strings)
        energyFile.close()

    except Exception as err:
        print(ERR_ENER_WRITE_FAIL)
        print('Energies to be written: (etot, efermi) =', energyPair)
        exit_with_error(err)
        
    return energyPair # A pair if it's useful for anything later

def get_energies(DIR_RELAXATION, DIR_ANALYSIS, writeOut=False):
    print('Fetching total energy calculations from calculations...')
    outcar = Outcar(DIR_RELAXATION + OUTCAR_NAME)
    mkdir(TOTAL_ENER_DIR_NAME, DIR_ANALYSIS)
    energy_pair = get_energy_analysis(checkPath(DIR_ANALYSIS + TOTAL_ENER_DIR_NAME), outcar, writeOut=writeOut) # energyPair = (total energy, fermi energy) tuple pair
    return energy_pair