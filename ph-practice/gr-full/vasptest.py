from pymatgen.io.vasp.outputs import Outcar
from ___constants_names import TOT_ENERGIES_NAME
from __directory_searchers import checkPath
from ___constants_misc import ERR_ENER_WRITE_FAIL

def get_energy_analysis(dirName, outcar):
    dirName = checkPath(dirName)
    try: 
        energyFile = open(dirName + TOT_ENERGIES_NAME, 'w')

        str1 = 'Total energy (eV): ' + str(outcar.final_energy)
        str2 = 'Fermi energy (eV): ' + str(outcar.efermi)

        strings = [str1, '\n', str2]

        energyFile.writelines(strings)
        energyFile.close()

        return (outcar.final_energy, outcar.efermi) # A pair if it's useful for anything later
    except Exception as err:
        print(ERR_ENER_WRITE_FAIL)
        print('Error:', err)

print(get_energy_analysis('.', Outcar('OUTCAR')))
