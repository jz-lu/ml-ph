import matplotlib.pyplot as plt
import matplotlib
import numpy as np
from ___constants_config import DEFAULT_ABS_MIN_ENERGY
from ___constants_names import DSAMPLE_ENERGIES_NAME, DSAMPLE_SPACINGS_NAME, DSAMPLE_ENERGIES_TXT, DSAMPLE_SPACINGS_TXT
from ___constants_output import DEFAULT_CONTOUR_LEVELS, HIGH_SYMMETRY_LABELS
from __directory_searchers import checkPath
from scipy.interpolate import make_interp_spline as interpolate_scatter

class DSamplingOutput:
    def __init__(self, out_dir, npts, special_pts=None, energies=None, spacings=None, ltype='hexagonal'):
        assert ltype == 'hexagonal', f"Only hexagonal lattices supported (for now), but got {ltype}"
        self.out_dir = checkPath(out_dir); self.npts = npts + 1 # add 1 for periodici boundary conditions
        minenergy = min(energies); self.energies = 1000*(np.array(energies)-minenergy)
        print(f"Adjusted energies (meV): {self.energies}")
        self.spacings = np.array(spacings)
        assert special_pts is None or len(special_pts) == 3, f"Must give list of 3 special points, but is {special_pts}"
        self.special_pts = np.array(special_pts); self.special_pts = np.append(self.special_pts, npts)
        self.energies = np.append(self.energies, self.energies[0])
        self.spacings = np.append(self.spacings, self.spacings[0]) # impose periodic boundary conditions
        print("Initialized DSamplingOutput object")
    def save_raw_data(self):
        np.save(self.out_dir + DSAMPLE_ENERGIES_NAME, self.energies)
        np.save(self.out_dir + DSAMPLE_SPACINGS_NAME, self.spacings)
        np.savetxt(self.out_dir + DSAMPLE_ENERGIES_TXT, self.energies)
        np.savetxt(self.out_dir + DSAMPLE_SPACINGS_TXT, self.spacings)
        print(f"Saved raw data over {self.npts-1} diagonally sampled points to {self.out_dir}")
    def __diag_plot(self, arr, plt_type='energy'):
        assert isinstance(plt_type, str); plt.clf(); fig, ax = plt.subplots()
        # ax.set_aspect('equal') # prevent axis stretching
        ax.set_title(f"{'Energies' if plt_type == 'energy' else 'Interlayer spacing'} along diagonal")
        x = np.linspace(0, 1, self.npts); y = self.energies if plt_type == 'energy' else self.spacings
        if self.special_pts is not None:
            print(f"Adding high-symmetry tick labels {HIGH_SYMMETRY_LABELS} at {x[self.special_pts]}")
            plt.xticks(ticks=x[self.special_pts], labels=HIGH_SYMMETRY_LABELS)
        else:
            plt.tick_params(
                axis='x',          # changes apply to the x-axis
                which='both',      # both major and minor ticks are affected
                bottom=False,      # ticks along the bottom edge are off
                top=False,         # ticks along the top edge are off
                labelbottom=False) # labels along the bottom edge are off
        ax.set_ylabel(r"$E_{tot} (meV)$" if plt_type == 'energy' else 'Interlayer spacing (unitless)')
        ax.scatter(x, y); fig.savefig(self.out_dir + f"diag_{plt_type}_scatter.png")
        line = ax.plot(x, y, c='k')
        fig.savefig(self.out_dir + f"diag_{plt_type}_jagged.png")
        line.pop(0).remove()
        interpol = interpolate_scatter(x, y)
        xp = np.linspace(0, 1, 301); yp = interpol(xp)
        ax.plot(xp, yp, c='k')
        fig.savefig(self.out_dir + f"diag_{plt_type}_smooth.png")
    def plot_energies(self):
        self.__diag_plot(self.energies, plt_type='energy')
    def plot_spacings(self):
        self.__diag_plot(self.spacings, plt_type='z')
    def output_all_analysis(self):
        self.save_raw_data(); self.plot_energies(); self.plot_spacings()

class ConfigOutput:
    # plot list is a list of tuples (b, z, e) = (shift, z-spacing, energy).
    # cob is the change-of-basis matrix to get fom lattice basis (which b is in) to Cartesian to plot.
    def __init__(self, out_dir, plot_list, cob_matrix, abs_min_energy=None):
        print("Initalizing ConfigOutput object.")
        # plot_list: list of (b, z, e) points
        self.__plot_list = plot_list
        self.__zspacings = np.array([i[1] for i in plot_list])
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

        print(f'Saving raw data to {self.__out_dir}')
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
        ax.set_title(r"Relaxed interlayer spacing (direct coordinates)")
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