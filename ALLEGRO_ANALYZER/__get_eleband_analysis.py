from pymatgen.io.vasp.outputs import Vasprun
from pymatgen import Spin
from pymatgen.electronic_structure.plotter import BSPlotter, BSDOSPlotter, DosPlotter
from ___constants_names import VASP_RUN_XML_NAME, PLOT_FILE_FORMAT, ELEBAND_RAW_DATA_NAME
from __dirModifications import checkPath
from __query_inputs import getInputFormula

def get_eleband_analysis(inputDir, outputDir, poscarObj, extract_raw_data=True, extract_plot=True):
    print('Starting electronic band structure parsing of vasp calculations...')
    inputDir = checkPath(inputDir)
    outputDir = checkPath(outputDir)

    # Get the full run data and extract band info
    run = Vasprun(inputDir + VASP_RUN_XML_NAME, parse_projected_eigen=True)
    structure = run.get_band_structure()
    band = BSPlotter(structure)

    if extract_raw_data:
        f = open(outputDir + ELEBAND_RAW_DATA_NAME, 'w')
        f.write(str(band.bs_plot_data(zero_to_efermi=True)))
        f.close()
        print('Electronic band structure raw data written to %s'%(outputDir))

    if extract_plot:
        name = getInputFormula(poscarObj) # Label
        band.save_plot(outputDir + '{}_eleband.'.format(name) + PLOT_FILE_FORMAT, img_format=PLOT_FILE_FORMAT)
        print('Electronic band structure plot written to %s'%(outputDir))

    return band

