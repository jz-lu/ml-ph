from pymatgen.electronic_structure.plotter import BSDOSPlotter
import matplotlib

from ____exit_with_error import exit_with_error

from __directory_searchers import checkPath

def get_elecombined_analysis(dirName, dosObj, bandObj):
    try:
        dirName = checkPath(dirName)
        plot = BSDOSPlotter.get_plot(bandObj, dosObj)

        plot.savefig(dirName)
        return plot
    except Exception as err:
        exit_with_error('Error in producing combined ELEDOS/ELEBAND plot: ' + err)