import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import numpy.linalg as LA
from ___constants_config import DEFAULT_ABS_MIN_ENERGY
from ___constants_names import (
    DSAMPLE_ENERGIES_NAME, DSAMPLE_SPACINGS_NAME, 
    DSAMPLE_ENERGIES_TXT, DSAMPLE_SPACINGS_TXT, 
    DSAMPLE_FORCES_NAME, 
    FGSFE_COEFF_NAME, FGSFE_COEFF_TXT, 
    FGSFE_SCORE_NAME, FGSFE_PLOT_NAME
)
from ___constants_output import DEFAULT_CONTOUR_LEVELS, HIGH_SYMMETRY_LABELS
from __directory_searchers import checkPath
from scipy.interpolate import make_interp_spline as interpolate_scatter
from pymatgen.io.vasp.inputs import Poscar
from sklearn.linear_model import LinearRegression
from math import pi, sqrt
from ___helpers_parsing import warn
import phonopy, copy, os

class DSamplingOutput:
    def __init__(self, out_dir, npts, special_pts=None, energies=None, spacings=None, ph_list=None, ltype='hexagonal'):
        assert ltype == 'hexagonal', f"Only hexagonal lattices supported (for now), but got {ltype}"
        print("Initializing DSamplingOutput object")
        self.out_dir = checkPath(out_dir); self.npts = npts + 1 # add 1 for periodici boundary conditions
        self.spacings = None; self.energies = None; self.force_matrices = None
        if energies is not None:
            minenergy = min(energies); self.energies = 1000*(np.array(energies)-minenergy)
            print(f"Adjusted energies (meV): {self.energies}")
            self.energies = np.append(self.energies, self.energies[0])
            print("Initialized energies")
        if spacings is not None:
            self.spacings = np.array(spacings)
            self.spacings = np.append(self.spacings, self.spacings[0]) # impose periodic boundary conditions
            print("Initialized interlayer spacings")
        assert special_pts is None or len(special_pts) == 3, f"Must give list of 3 special points, but is {special_pts}"
        if special_pts is not None:
            self.special_pts = np.array(special_pts); self.special_pts = np.append(self.special_pts, npts)
            print("Initialized special points")
        if ph_list is not None:
            self.force_matrices = [ph.get_dynamical_matrix_at_q([0,0,0]) for ph in ph_list]
            print("Initialized phonopy object list")

    def save_raw_data(self):
        if self.energies is not None:
            np.save(self.out_dir + DSAMPLE_ENERGIES_NAME, self.energies)
            np.savetxt(self.out_dir + DSAMPLE_ENERGIES_TXT, self.energies)
        if self.spacings is not None:
            np.save(self.out_dir + DSAMPLE_SPACINGS_NAME, self.spacings)
            np.savetxt(self.out_dir + DSAMPLE_SPACINGS_TXT, self.spacings)
        if self.force_matrices is not None:
            np.savez(self.out_dir + DSAMPLE_FORCES_NAME, *self.force_matrices)
        print(f"Saved raw data over {self.npts-1} diagonally sampled points to {self.out_dir}")

    def __diag_plot(self, arr, plt_type='energy'):
        assert plt_type in ['energy', 'z', 'forces']; plt.clf(); fig, ax = plt.subplots()
        title = 'Energies'; y_lab = r"$E_{tot} (meV)$"; y = self.energies
        if plt_type == 'z':
            title = 'Interlayer spacing'; y_lab = 'Interlayer spacing (unitless)'
        elif plt_type == 'forces':
            pass
        ax.set_title(f"{title} along diagonal")
        x = np.linspace(0, 1, self.npts); y = self.energies if plt_type == 'energy' else self.spacings
        print(f"Now plotting TYPE={plt_type}")
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
        ax.set_ylabel(y_lab)
        ax.scatter(x, y); fig.savefig(self.out_dir + f"diag_{plt_type}_scatter.png")
        interpol = interpolate_scatter(x, y)
        xp = np.linspace(0, 1, 301); yp = interpol(xp)
        ax.plot(xp, yp, c='k')
        fig.savefig(self.out_dir + f"diag_{plt_type}_smooth.png")
    
    def plot_energies(self):
        assert self.energies is not None, "Energy data was not provided to analyzer"
        self.__diag_plot(self.energies, plt_type='energy')
    
    def plot_spacings(self):
        assert self.spacings is not None, "Interlayer spacings data was not provided to analyzer"
        self.__diag_plot(self.spacings, plt_type='z')
    
    def output_all_analysis(self):
        self.save_raw_data(); self.plot_energies(); self.plot_spacings()
    
    def plot_forces(self, atomic_idx_pairs, layer_idx_pairs, cart_idxs, poscar : Poscar):
        cols = list(mcolors.TABLEAU_COLORS.keys()); ncol = len(cols)
        assert len(atomic_idx_pairs) == len(layer_idx_pairs) == len(cart_idxs), f"[Set 1] Number of atomic, layer, and Cartesian indices must be the same, but is {len(atomic_idx_pairs)}, {len(layer_idx_pairs)}, and {len(cart_idxs)}"
        atomic_sites = list(map(lambda x: x.species.elements[0].symbol, poscar.structure.sites)); n_at = len(atomic_sites)
        cart_letters = copy.deepcopy(cart_idxs)
        for i, idx in enumerate(cart_idxs):
                if idx == 'x':
                    cart_idxs[i] = 0
                elif idx == 'y':
                    cart_idxs[i] = 1
                elif idx == 'z':
                    cart_idxs[i] = 2
                else:
                    assert False, f"Invalid Cartesian element {idx} at position {i}"
        print(f"Transformed Cartesian pairs to {cart_idxs}")
        plt.clf(); fig, ax = plt.subplots(); x = np.linspace(0, 1, self.npts)
        ax.set_title(f"Forces along diagonal")
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
        ax.set_ylabel(r"Force (ev/$\AA$)")
        for i, (ats, ls, c, cl) in enumerate(zip(atomic_idx_pairs, layer_idx_pairs, cart_idxs, cart_letters)):
            ats = np.array(ats); ls = np.array(ls); c = np.array(c)
            assert len(ats) == len(ls) == len(c) == 2
            assert ls[0] in [1, 2] and ls[1] in [1,2]; assert 0 <= ats[0] < n_at and 0 <= ats[1] < n_at
            idxs = 3*ats + c; y = np.array([f[idxs[0], idxs[1]] for f in self.force_matrices])
            interpol = interpolate_scatter(x, y); xp = np.linspace(0, 1, 301); yp = interpol(xp)
            lab = '%s%d-%s%d(%x)'%(atomic_sites[ats[0]], ls[0], atomic_sites[ats[1]], ls[1], cl)
            ax.scatter(x, y, c=cols[i%ncol]); ax.plot(xp, yp, c=cols[i%ncol], label=lab)
        if len(cart_letters) > ncol:
            warn("Warning: there are more force lines than colors available. Some lines will be ambiguous.")
        ax.legend()
        fig.savefig(self.out_dir + f"diag_forces.png")

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
        self.__direct_shifts = [np.array(i[0][:2]) for i in plot_list]
        self.b1shifts = np.array([i[0] for i in self.__direct_shifts]) # direct coordinates
        self.b2shifts = np.array([i[1] for i in self.__direct_shifts])
        self.__shifts = [np.dot(cob_matrix, i) for i in self.__direct_shifts] # Cartesian coordinates
        self.xshifts = np.array([i[0] for i in self.__shifts])
        self.yshifts = np.array([i[1] for i in self.__shifts])
        print("DIRECT SHIFTS:", self.__direct_shifts)
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
    
# Use basis linear regression to fit GSFE to leading 6 fourier terms
# See Eq(4), Carr 2018
# Bugs: may not work if unit cell vectors is 60 deg instead of 120
class FourierGSFE:
    def __init__(self, energies, b_matrix, lattice_matrix, ltype='hexagonal', sampling_type='grid'):
        self.stype = sampling_type; assert isinstance(self.stype, str)
        self.__A = lattice_matrix
        print("Initializing FourierGSFE object")
        assert ltype == 'hexagonal', "Non-hexagonal lattices not supported"
        assert len(energies) == len(b_matrix), f"Must have same number of configurations as energies, but got {len(b_matrix)} vs. {len(energies)}"
        b_matrix = b_matrix[:,:2]
        print(f"Configurations:\n {b_matrix}")
        minenergy = min(energies); self.GSFE = 1000*(np.array(energies)-minenergy)
        print(f"Adjusted energies (meV): {self.GSFE}")
        self.nb = len(b_matrix); self.b_matrix = b_matrix
        self.__fitted = False; self.coeffs = None; self.reg = None
        self.X = self.__b_to_fourier_basis(b_matrix)
    def __b_to_vw(self, b_mat):
        M = 2 * pi * np.array([[1, -1/sqrt(3)], [0, 2/sqrt(3)]]) / LA.norm(self.__A[0]) # for hexagonal lattices of 120 lattice angle only
        return (M @ b_mat.T).T
    def __vw_to_fourier_basis(self, vw):
        X = np.ones((self.nb, 5)); v = vw[:,0]; w = vw[:,1] # col 0 is bias
        X[:,0] = np.cos(v) + np.cos(w) + np.cos(v + w)
        X[:,1] = np.cos(v + 2*w) + np.cos(v - w) + np.cos(2*v + w) 
        X[:,2] = np.cos(2*v) + np.cos(2*w) + np.cos(2*v + 2*w) 
        X[:,3] = np.sin(v) + np.sin(w) - np.sin(v + w)
        X[:,4] = np.sin(2*v + 2*w) - np.sin(2*v) - np.sin(2*w)
        return X
    def __b_to_fourier_basis(self, b_mat):
        return self.__vw_to_fourier_basis(self.__b_to_vw(b_mat))
    def __fit_coeff(self):
        assert not self.__fitted, f"Fitting already done"
        print("Fitting energies to Fourier series...")
        self.reg = LinearRegression().fit(self.X, self.GSFE)
        self.coeffs = np.append(self.reg.intercept_, self.reg.coef_)
        self.__fitted = True
        return self.coeffs
    def __ensure_fitted(self):
        if not self.__fitted:
            self.__fit_coeff()
        assert self.coeffs is not None and self.reg is not None, f"Fitting failed"
    def save_raw_data(self, outdir):
        print("Saving coefficents and score to file...")
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
        self.__ensure_fitted()
        np.save(outdir + FGSFE_COEFF_NAME, self.coeffs)
        np.save(outdir + FGSFE_COEFF_TXT, self.coeffs)
        np.save(outdir + FGSFE_SCORE_NAME, self.get_score())
    def get_coeffs(self):
        self.__ensure_fitted(); return self.coeffs
    def get_score(self):
        self.__ensure_fitted(); return self.reg.score(self.X, self.GSFE)
    def predict(self, b_matrix):
        self.__ensure_fitted(); return self.reg.predict(self.__b_to_fourier_basis(b_matrix))
    def plot_pred_vs_actual(self, outdir, outname=FGSFE_PLOT_NAME):
        print("Plotting predicted vs. actual...")
        assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
        outdir = checkPath(outdir); plt.clf(); fig, ax = plt.subplots()
        x = self.predict(self.b_matrix); y = self.GSFE
        maxmax = max(max(x), max(y)); minmin = min(min(x), min(y))
        ax.set_xlabel("Predicted"); ax.set_ylabel("Actual")
        ax.scatter(x, y, c='k')
        ax.set_title(f"GSFE Fourier fitting on {self.stype} sampling")
        ax.plot([minmin, maxmax], [minmin, maxmax], c='royalblue') # y=x line
        print(f"Predicted:\n {x}\nActual: \n{y}")
        fig.savefig(outdir + outname)
    def output_all_analysis(self, outdir):
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
        self.__ensure_fitted()
        self.save_raw_data(outdir)
        self.plot_pred_vs_actual(outdir)
        print("All Fourier GSFE analysis written to file")

