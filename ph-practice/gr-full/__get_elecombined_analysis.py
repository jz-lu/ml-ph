from pymatgen.electronic_structure.plotter import BSDOSPlotter
import matplotlib

from __directory_searchers import checkPath

def get_elecombined_analysis(dirName, dosObj, bandObj):
    dirName = checkPath(dirName)
    plot = BSDOSPlotter.get_plot(bandObj, dosObj)

    plot.savefig(dirName)
    return plot