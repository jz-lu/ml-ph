import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from __directory_searchers import checkPath

NUM_LEVELS = 21 # Controls contour plot smooothness, ensure it is odd for the middle to be captured.

class DataOutput:
    # plot list is a list of tuples (b, z, e) = (shift, z-spacing, energy).
    # cob is the change-of-basis matrix to get fom lattice basis (which b is in) to Cartesian to plot.
    #   Note cob_matrix must be a numpy matrix.
    def __init__(self, out_dir, plot_list, cob_matrix):
        print("Initalizing DataOutput object.")
        self.__plot_list = plot_list
        self.__zspacings = np.array([plot_list[1] for i in plot_list])
        self.__energies = np.array([plot_list[2] for i in plot_list])
        self.__out_dir = checkPath(out_dir)
        print("Output to be stored in %s"%(out_dir))

        # Get the shifts in Cartesian coordinates via COB
        self.__shifts = [np.dot(cob_matrix, np.array(plot_list[0][:-1])) for i in plot_list] # Exclude z-coordinate since plot is 2D
        self.xshifts = np.array([i[0] for i in self.__shifts])
        self.yshifts = np.array([i[1] for i in self.__shifts])
        print("SHIFTS:", self.__shifts)
        print("xshifts:", self.xshifts)
        print("yshifts:", self.yshifts)
        print("Energies:", self.__energies)
        print("Interlayer spacings:", self.__zspacings)
    
    # Output raw data as a csv file.
    def save_raw_data(self):
        table = np.transpose(np.array([self.xshifts, self.yshifts, self.__zspacings, self.__energies]))
        out_file = self.__out_dir + "raw_data.csv"
        with open(out_file) as f:
            f.write(b"2D Shift vector (lattice basis), Relaxed interlayer spacing, Energy\n")
            np.savetxt(out_file, table, delimiter=",")
   
    # Smoothly interpolate toal energy at various shifts.
    def plot_e_vs_b(self):
        fig, ax = plt.subplots()
        cf = ax.tricontourf(self.xshifts, self.yshifts, self.__energies, levels=NUM_LEVELS, cmap="RdGy")
        fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$b_x$")
        ax.set_ylabel(r"$b_y$")
        ax.set_title(r"$E_{tot}(b)$")
        out_file = self.__out_dir + "energy_config_plot.png"
        fig.savefig(out_file)
    
    # Smoothly interpolate interlayer spacing at various shifts.
    def plot_z_vs_b(self):
        fig, ax = plt.subplots()
        cf = ax.tricontourf(self.xshifts, self.yshifts, self.__zspacings, levels=NUM_LEVELS, cmap="twilight_shifted")
        fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$b_x$")
        ax.set_ylabel(r"$b_y$")
        ax.set_title("Relaxed interlayer spacing (z)")
        out_file = self.__out_dir + "z_config_plot.png"
        fig.savefig(out_file)
    
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
    
    # Do all available functions.
    def output_all_analysis(self):
        self.save_raw_data()
        self.plot_e_vs_z()
        self.plot_e_vs_b()
        self.plot_z_vs_b()