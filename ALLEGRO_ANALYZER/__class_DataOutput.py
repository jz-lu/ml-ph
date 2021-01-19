import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from __directory_searchers import checkPath

class DataOutput:
    def __init__(self, out_dir, plot_list): # plot list is a list of tuples (b, z, e) = (shift, z-spacing, energy)
        self.__plot_list = plot_list
        self.__shifts = [plot_list[0] for i in plot_list]
        self.__zspacings = [plot_list[1] for i in plot_list]
        self.__energies = [plot_list[2] for i in plot_list]
        self.__out_dir = checkPath(out_dir)
    
    # Output raw data as a csv file.
    def save_raw_data(self):
        table = np.transpose(np.array([self.__shifts, self.__zspacings, self.__energies]))
        out_file = self.__out_dir + "raw_data.csv"
        with open(out_file) as f:
            f.write(b"Shift vector (lattice basis), Relaxed interlayer spacing, Energy\n")
            np.savetxt(out_file, table, delimiter=",")

    
    # Smoothly interpolate toal energy at various shifts.
    def plot_e_vs_b(self):
        # TODO
    
    # Smoothly interpolate interlayer spacing at various shifts.
    def plot_z_vs_b(self):
        # TODO
    
    # Evaluate energy vs. interlayer spacing over all of the shifts.
    def plot_e_vs_z(self):
        # Line plot
        fig, ax = plt.subplots()
        ax.plot(self.__zspacings, self.__energies)
        fig.savefig(self.__out_dir + "energy_vs_interlayer_spacing_line.png")
        
        # Scatter plot
        scat, ax2 = plt.subplots()
        ax2.scatter(self.__zspacings, self.__energies)
        scat.savefig(self.__out_dir + "energy_vs_interlayer_spacing_scatter.png")
        