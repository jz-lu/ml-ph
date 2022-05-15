from pymatgen.io.vasp.outputs import Vasprun
from pymatgen import Spin
from pymatgen.electronic_structure.plotter import DosPlotter

from ___constants_names import VASP_RUN_XML_NAME, ELEDOS_RAW_DATA_NAME, PLOT_FILE_FORMAT

from __directory_searchers import checkPath
from __query_inputs import getInputFormula

def get_eledos_analysis(inputDir, outputDir, poscarObj, extract_raw_data=True, extract_plot=True):
    print('Starting parsing of electronic DOS data from vasp run...')
    inputDir = checkPath(inputDir)
    outputDir = checkPath(outputDir)

    # Get the full run data and extract DOS
    run = Vasprun(inputDir + VASP_RUN_XML_NAME)
    dos = run.tdos

    # Just name by the number of atoms
    name = getInputFormula(poscarObj)

    # Get the DOS info
    dosplot = DosPlotter(stack=True)
    dosplot.add_dos("Total DOS for {}".format(str(name)), dos)

    if extract_raw_data:
        f = open(outputDir + ELEDOS_RAW_DATA_NAME, 'w')
        f.write(str(dosplot.get_dos_dict()))
        f.close()
        print('Electronic DOS raw data written to %s'%(outputDir))

    if extract_plot:
        dosplot.save_plot(outputDir + '{}_eledos_wide.'.format(name) + PLOT_FILE_FORMAT, PLOT_FILE_FORMAT)
        dosplot.save_plot(outputDir + '{}_eledos_narrow.'.format(name) + PLOT_FILE_FORMAT, PLOT_FILE_FORMAT, [-5, 5])
        print('Electronic DOS plot written to %s'%(outputDir))
    
    return dos



