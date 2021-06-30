import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from ___constants_config import DEFAULT_ABS_MIN_ENERGY
from ___constants_output import DEFAULT_CONTOUR_LEVELS
from __directory_searchers import checkPath

class DataOutput:
    # plot list is a list of tuples (b, z, e) = (shift, z-spacing, energy).
    # cob is the change-of-basis matrix to get fom lattice basis (which b is in) to Cartesian to plot.
    def __init__(self, out_dir, plot_list, cob_matrix, abs_min_energy=None):
        print("Initalizing DataOutput object.")
        # plot_list: list of (b, z, e) points
        self.__plot_list = plot_list
        minz = min([i[1] for i in plot_list])
        self.__zspacings = np.array([(i[1]-minz)*1000 for i in plot_list])
        if not abs_min_energy:
            abs_min_energy = min([i[2] for i in plot_list])
        print(f"Shifting by minimum energy {abs_min_energy} eV")
        self.__energies = np.array([(i[2]-abs_min_energy)*1000 for i in plot_list])
        self.__out_dir = checkPath(out_dir)
        print("Output to be stored in %s"%(out_dir))

        # Get the shifts in Cartesian coordinates via COB
        self.__direct_shifts = [np.array(i[0][:-1]) for i in plot_list]
        self.b1shifts = np.array([i[0] for i in self.__direct_shifts]) # direct coordinates
        self.b2shifts = np.array([i[1] for i in self.__direct_shifts])
        self.__shifts = [np.dot(cob_matrix, i) for i in self.__direct_shifts] # Cartesian coordinates
        self.xshifts = np.array([i[0] for i in self.__shifts])
        self.yshifts = np.array([i[1] for i in self.__shifts])
        print("DIRECT SHIFTS:", [np.array(i[0][:-1]) for i in plot_list])
        print("CARTESIAN SHIFTS:", self.__shifts)
        print("xshifts:", self.xshifts)
        print("yshifts:", self.yshifts)
        print("Energies (meV):", self.__energies)
        print("Interlayer spacings:", self.__zspacings)
        np.save(self.__out_dir + 'bze', self.__plot_list)

        print('[DEBUG] RAW DATA OUTPUTTING!')
        with open(self.__out_dir + "shifts.txt", 'w+') as f1:
            f1.write(str(self.__plot_list))
        with open(self.__out_dir + "xshifts.txt", 'w+') as f1:
            f1.write(str(self.xshifts))
        with open(self.__out_dir + "yshifts.txt", 'w+') as f1:
            f1.write(str(self.yshifts))
        with open(self.__out_dir + "zspacings.txt", 'w+') as f1:
            f1.write(str(self.__zspacings))
        with open(self.__out_dir + "e.txt", 'w+') as f1:
            f1.write(str(self.__energies))

    # Output raw data as a csv file.
    def save_raw_data(self):
        table = np.transpose(np.array([self.xshifts, self.yshifts, self.__zspacings, self.__energies]))
        out_file = self.__out_dir + "raw_data.csv"
        np.savetxt(out_file, table, delimiter=",", 
            header="bx, by, relaxed z-spacing, energy\n")
   
    # Smoothly interpolate toal energy at various shifts.
    def plot_e_vs_b(self, levels=DEFAULT_CONTOUR_LEVELS):
        plt.clf(); fig, ax = plt.subplots()
        cf = ax.tricontourf(self.xshifts, self.yshifts, self.__energies, 
                            levels=levels, cmap="RdGy")
        fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$b_x$")
        ax.set_ylabel(r"$b_y$")
        ax.set_title(r"$E_{tot}(\mathbf{b}) (meV)$")
        out_file = self.__out_dir + f"energy_config_plot_cart_{levels}"
        ax.set_aspect('equal') # prevent axis stretching
        fig.savefig(out_file + "_eqasp.png")
        plt.clf(); fig, ax = plt.subplots()
        cf = ax.tricontourf(self.b1shifts, self.b2shifts, self.__energies, 
                            levels=levels, cmap="RdGy"); fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$a_1$")
        ax.set_ylabel(r"$a_2$")
        ax.set_title(r"$E_{tot}(\mathbf{b}=b_1 \mathbf{a}_1 + b_2 \mathbf{a}_2) (meV)$")
        out_file = self.__out_dir + f"energy_config_plot_direct_{levels}"
        fig.savefig(out_file + "_eqasp.png")
    
    # Smoothly interpolate interlayer spacing at various shifts.
    def plot_z_vs_b(self, levels=DEFAULT_CONTOUR_LEVELS,):
        fig, ax = plt.subplots()
        cf = ax.tricontourf(self.xshifts, self.yshifts, self.__zspacings, 
                            levels=levels, cmap="RdGy") # or try cmap="twilight_shifted" for a fun twist
        fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$b_x$")
        ax.set_ylabel(r"$b_y$")
        ax.set_title("Relaxed interlayer spacing (normed units)")
        out_file = self.__out_dir + "z_config_plot_cart"
        ax.set_aspect('equal') # prevent axis stretching
        fig.savefig(out_file + "_eqasp.png")
    
    # Evaluate energy vs. interlayer spacing over all of the shifts.
    def plot_e_vs_z(self):
        scat, ax2 = plt.subplots()
        ax2.scatter(self.__zspacings, self.__energies)
        scat.savefig(self.__out_dir + "energy_vs_interlayer_spacing_scatter.png")
    
    # Do all available functions.
    def output_all_analysis(self, levels=DEFAULT_ABS_MIN_ENERGY):
        self.save_raw_data()
        self.plot_e_vs_b(levels=levels)
        self.plot_z_vs_b(levels=levels)
        self.plot_e_vs_z()