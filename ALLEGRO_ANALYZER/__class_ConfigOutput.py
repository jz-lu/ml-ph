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
    def __init__(self, out_dir, npts, name, special_pts=None, energies=None, spacings=None, scaled=False, dump=True, 
                 ph_list=None, ltype='hexagonal'):
        assert ltype == 'hexagonal', f"Only hexagonal lattices supported (for now), but got {ltype}"
        print("Initializing DSamplingOutput object")
        assert npts % 3 == 0, 'Must choose a number of points as a multiple of 3 to respect symmetry'
        self.out_dir = checkPath(out_dir); self.npts = npts + 1 # add 1 for periodici boundary conditions
        self.spacings = None; self.energies = None; self.force_matrices = None
        self.name = name
        if energies is not None:
            minenergy = min(energies); self.energies = (1 if scaled else 1000)*(np.array(energies)-minenergy)
            if dump:
                print(f"Adjusted energies (meV): {self.energies}")
            self.energies = np.append(self.energies, self.energies[0])
            assert len(self.energies) == self.npts, f"{len(self.energies)} and {self.npts} points"
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

    def __diag_plot(self, arr, plt_type='energy', pfx='', tsfx='', 
                    interp=True, scat=True, line=False, addendum=None):
        assert interp or scat or line, "Must specify at least 1 type of plot"
        assert not (interp and line), "Must choose at most one of interpolation and line plotting"
        assert plt_type in ['energy', 'z', 'forces']; plt.clf(); fig, ax = plt.subplots()
        title = 'GSFE'; y_lab = r"$E_{tot} (meV)$"; y = self.energies
        if plt_type == 'z':
            title = 'Interlayer spacing'; y_lab = 'Interlayer spacing (unitless)'
        title = title + (f' of {tsfx}' if tsfx != '' else '')
        print(f"Title: {title} (sfx={tsfx})")
        ax.set_title(f"{title} along diagonal")
        x = np.linspace(0, 1, self.npts); y = self.energies if plt_type == 'energy' else self.spacings
        print(f"Now plotting TYPE={plt_type} of {len(x)} points")
        if self.special_pts is not None:
            print(f"Adding high-symmetry tick labels {HIGH_SYMMETRY_LABELS} at {x[self.special_pts]}")
            plt.xticks(ticks=x[self.special_pts], labels=HIGH_SYMMETRY_LABELS)
        else:
            plt.tick_params(
                axis='x',           # changes apply to the x-axis
                which='both',       # both major and minor ticks are affected
                bottom=False,       # ticks along the bottom edge are off
                top=False,          # ticks along the top edge are off
                labelbottom=False)  # labels along the bottom edge are off
        ax.set_ylabel(y_lab)
        if scat:
            print(f"Scattering x,y with {len(x)} points")
            ax.scatter(x, y, c='darkslategrey')
        if addendum is not None:
            print(f"Scattering addendum {len(addendum[0])} points")
            ax.scatter(addendum[0], addendum[1], c='darkslategrey')
        if interp:
            interpol = interpolate_scatter(x, y)
            xp = np.linspace(0, 1, 301); yp = interpol(xp)
            ax.plot(xp, yp, c='royalblue')
        if line:
            print(f"Lining x,y with {len(x)} points")
            ax.plot(x, y, c='royalblue')
        fig.savefig(self.out_dir + pfx + f"diag_{plt_type}.png")
    
    def plot_energies(self, pfx='', tsfx='', interp=True, scat=True, line=False, addendum=None):
        assert self.energies is not None, "Energy data was not provided to analyzer"
        self.__diag_plot(self.energies, plt_type='energy', pfx=pfx, 
                         tsfx=tsfx, interp=interp, scat=scat, line=line, addendum=addendum)
    
    def plot_spacings(self, pfx='', interp=True, scat=True, line=False):
        assert self.spacings is not None, "Interlayer spacings data was not provided to analyzer"
        self.__diag_plot(self.spacings, plt_type='z', pfx=pfx, interp=interp, scat=scat, line=line)
    
    def output_all_analysis(self):
        self.save_raw_data(); self.plot_energies(); self.plot_spacings()
    
    def plot_forces(self, atomic_idx_pairs, layer_idx_pairs, cart_idx_pairs, poscar : Poscar, pfx=''):
        print(f"Atomic pairs: {atomic_idx_pairs}, layer pairs: {layer_idx_pairs}, cart pairs: {cart_idx_pairs}")
        cols = list(mcolors.TABLEAU_COLORS.keys()); ncol = len(cols)
        assert self.force_matrices is not None, f"Must give valid list of phonopy objects to plot forces"
        assert len(atomic_idx_pairs) == len(layer_idx_pairs) == len(cart_idx_pairs), f"Number of atomic, layer, and Cartesian indices must be the same, but is {len(atomic_idx_pairs)}, {len(layer_idx_pairs)}, and {len(cart_idx_pairs)}"
        atomic_sites = list(map(lambda x: x.species.elements[0].symbol, poscar.structure.sites)); n_at = len(atomic_sites)
        nplt = len(atomic_idx_pairs); assert int(sqrt(nplt))**2 == nplt, f"Number to plot {nplt} must be a perfect square to plot properly"
        cart_letter_pairs = [['']*cart_idx_pairs.shape[1] for _ in range(cart_idx_pairs.shape[0])]
        for i, pair in enumerate(cart_idx_pairs):
            for j, idx in enumerate(pair):
                if idx == 0:
                    cart_letter_pairs[i][j] = 'x'
                elif idx == 1:
                    cart_letter_pairs[i][j] = 'y'
                elif idx == 2:
                    cart_letter_pairs[i][j] = 'z'
                else:
                    assert False, f"Invalid Cartesian element {idx} at position {i}"
        print(f"Transformed Cartesian pairs to letter form: \n{cart_letter_pairs}")
        rows = int(sqrt(nplt))
        plt.clf(); fig = plt.figure(); fig.subplots(nrows=rows, ncols=rows) 
        x = np.linspace(0, 1, self.npts)
        plt.suptitle(f"Mass-scaled forces along diagonal")
        fig.text(0, 0.5, r"Force constants (eV$\cdot amu/\AA^2$)", va='center', rotation='vertical')
        for i, (ats, ls, cs, cl) in enumerate(zip(atomic_idx_pairs, layer_idx_pairs, cart_idx_pairs, cart_letter_pairs)):
            ats = np.array(ats); ls = np.array(ls); cs = np.array(cs)
            assert len(ats) == len(ls) == len(cs) == 2
            assert ls[0] in [1, 2] and ls[1] in [1,2]; assert 0 <= ats[0] < n_at and 0 <= ats[1] < n_at
            idxs = 3*ats + cs; y = np.array([f[idxs[0], idxs[1]] for f in self.force_matrices])
            y = np.append(y, y[0]) # impose periodic boundary conditions
            # interpol = interpolate_scatter(x, y); xp = np.linspace(0, 1, 301); yp = interpol(xp)
            lab = r'$%s^{(%d)}_%s \sim %s^{(%d)}_%s$'%(atomic_sites[ats[0]], ls[0], cl[0], atomic_sites[ats[1]], ls[1], cl[1])
            fig.axes[i].scatter(x, y, c=cols[i%ncol])
            fig.axes[i].text(0.7, 1.05, lab, transform=fig.axes[i].transAxes, size=10, weight='ultralight')
            fig.axes[i].ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            if self.special_pts is not None:
                plt.setp(fig.axes, xticks=x[self.special_pts], xticklabels=HIGH_SYMMETRY_LABELS)
            else:
                fig.axes[i].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
        plt.tight_layout()
        fig.savefig(self.out_dir + f"diag_forces{'' if pfx == '' else '_'}{pfx}.png")

class ConfigOutput:
    # plot list is a list of tuples (b, z, e) = (shift, z-spacing, energy).
    # cob is the change-of-basis matrix to get fom lattice basis (which b is in) to Cartesian to plot.
    def __init__(self, out_dir, plot_list, cob_matrix, name, 
                 ph_list=None, abs_min_energy=None, scaled=False, dump=False):
        self.name = name
        print(f"Initalizing ConfigOutput object on {self.name}.")
        self.cob = cob_matrix; self.ph_list = ph_list
        # plot_list: list of (b, z, e) points
        self.__plot_list = plot_list
        self.__zspacings = np.array([i[1] for i in plot_list])
        if not abs_min_energy:
            abs_min_energy = min([i[2] for i in plot_list])
        print(f"Shifting by minimum energy {abs_min_energy} eV")
        self.__energies = np.array([(i[2]-abs_min_energy)*(1 if scaled else 1000) for i in plot_list])
        self.__out_dir = checkPath(out_dir)
        print("Output to be stored in %s"%(out_dir))

        # Get the shifts in Cartesian coordinates via COB
        self.__direct_shifts = [np.array(i[0][:2]) for i in plot_list]
        self.b1shifts = np.array([i[0] for i in self.__direct_shifts]) # direct coordinates
        self.b2shifts = np.array([i[1] for i in self.__direct_shifts])
        self.__shifts = [np.dot(cob_matrix, i) for i in self.__direct_shifts] # Cartesian coordinates
        self.xshifts = np.array([i[0] for i in self.__shifts])
        self.yshifts = np.array([i[1] for i in self.__shifts])
        if dump:
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
        np.save(self.__out_dir + 'e.npy', self.__energies)

    def plot_diag_cut(self, energies, bprime, pfx='cut'):
        """
        Function assumes lattice angle is 60 deg.
        """
        bdir = np.array(self.__direct_shifts)
        diags = [np.isclose(i[0], i[1]) for i in bdir]
        print(f"Diags:\n{bdir[diags]}")
        e = self.__energies[diags]; e = np.append(e, e[0])
        addendum = [np.linspace(0,1,sum(diags)+1), e]
        diags = [np.isclose(i[0], i[1]) for i in bprime]
        e = energies[diags]
        ndiag = sum(diags)
        do = DSamplingOutput(self.__out_dir, ndiag, self.name, 
                        special_pts=[0, ndiag//3, 2*ndiag//3], energies=e, scaled=True)
        do.plot_energies(interp=False, line=True, addendum=addendum, pfx=pfx, scat=False, tsfx=self.name)
        
   
    # Smoothly interpolate toal energy at various shifts.
    def plot_e_vs_b(self, levels=DEFAULT_CONTOUR_LEVELS, pfx='', tpfx='', energies=None, b=None):
        need_direct = True
        if energies is None:
            energies = self.__energies
        else:
            assert b is not None
            need_direct = False
        if b is None:
            b = np.array(self.__shifts)
        else:
            assert energies is not None
            b = (self.cob @ np.array(b)[:,:2].T).T
        bx = b[:,0]; by = b[:,1]
        bdir = np.array(self.__direct_shifts)
        plt.clf(); fig, ax = plt.subplots()
        cf = ax.tricontourf(bx, by, energies, levels=levels, cmap="RdGy")
        fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$b_x$")
        ax.set_ylabel(r"$b_y$")
        ax.set_title(tpfx + ('' if tpfx=='' else ' ') + r"GSFE$(\mathbf{b})$ (meV) of %s"%self.name)
        out_file = self.__out_dir + pfx + f"energy_config_plot_{levels}.png"
        ax.set_aspect('equal') # prevent axis stretching
        fig.savefig(out_file)
        if need_direct:
            plt.clf(); fig, ax = plt.subplots()
            cf = ax.tricontourf(bdir[:,0], bdir[:,1], energies, levels=levels, cmap="RdGy"); fig.colorbar(cf, ax=ax)
            ax.set_xlabel(r"$a_1$")
            ax.set_ylabel(r"$a_2$")
            ax.set_title(r"$E_{tot}(\mathbf{b}=b_1 \mathbf{a}_1 + b_2 \mathbf{a}_2)$ (meV) of %s"%self.name)
            out_file = self.__out_dir + pfx + f"energy_config_plot_direct_{levels}.png"
            fig.savefig(out_file)
    
    # Smoothly interpolate interlayer spacing at various shifts.
    def plot_z_vs_b(self, levels=DEFAULT_CONTOUR_LEVELS,):
        fig, ax = plt.subplots()
        cf = ax.tricontourf(self.xshifts, self.yshifts, self.__zspacings, 
                            levels=levels, cmap="RdGy") # or try cmap="twilight_shifted" for a fun twist
        fig.colorbar(cf, ax=ax)
        ax.set_xlabel(r"$b_x$")
        ax.set_ylabel(r"$b_y$")
        ax.set_title(f"Relaxed interlayer spacing (direct) of {self.name}")
        out_file = self.__out_dir + "z_config_plot_cart"
        ax.set_aspect('equal') # prevent axis stretching
        fig.savefig(out_file + "_eqasp.png")
    
    # Evaluate energy vs. interlayer spacing over all of the shifts.
    def plot_e_vs_z(self):
        scat, ax2 = plt.subplots()
        ax2.scatter(self.__zspacings, self.__energies)
        scat.savefig(self.__out_dir + "energy_vs_interlayer_spacing_scatter.png")
    
    def plot_forces(self, atomic_idx_pairs, layer_idx_pairs, 
                    cart_idx_pairs, poscar : Poscar, 
                    levels=DEFAULT_CONTOUR_LEVELS, pfx=''):
        assert self.ph_list is not None, f"Must give list of phonopy objects to plot forces"
        self.force_matrices = [np.real(ph.get_dynamical_matrix_at_q([0,0,0])) for ph in self.ph_list]
        # cols = list(mcolors.TABLEAU_COLORS.keys()); ncol = len(cols)
        assert len(atomic_idx_pairs) == len(layer_idx_pairs) == len(cart_idx_pairs), f"Number of atomic, layer, and Cartesian indices must be the same, but is {len(atomic_idx_pairs)}, {len(layer_idx_pairs)}, and {len(cart_idx_pairs)}"
        atomic_sites = list(map(lambda x: x.species.elements[0].symbol, poscar.structure.sites)); n_at = len(atomic_sites)
        nplt = len(atomic_idx_pairs); assert int(sqrt(nplt))**2 == nplt, f"Number to plot {nplt} must be a perfect square to plot properly"
        cart_letter_pairs = [['']*cart_idx_pairs.shape[1] for _ in range(cart_idx_pairs.shape[0])]
        for i, pair in enumerate(cart_idx_pairs):
            for j, idx in enumerate(pair):
                if idx == 0:
                    cart_letter_pairs[i][j] = 'x'
                elif idx == 1:
                    cart_letter_pairs[i][j] = 'y'
                elif idx == 2:
                    cart_letter_pairs[i][j] = 'z'
                else:
                    assert False, f"Invalid Cartesian element {idx} at position {i}"
        print(f"Transformed Cartesian pairs to letter form: \n{cart_letter_pairs}")
        rows = int(sqrt(nplt))
        b = self.__shifts; bx = [i[0] for i in b]; by = [i[1] for i in b]
        plt.clf(); fig = plt.figure(); fig.subplots(nrows=rows, ncols=rows, sharex='col', sharey='row') 
        plt.suptitle(f"Forces over configurations")
        fig.text(0, 0.5, r"Force constants (eV/$\AA^2$)", va='center', rotation='vertical')
        for i, (ats, ls, cs, cl) in enumerate(zip(atomic_idx_pairs, layer_idx_pairs, cart_idx_pairs, cart_letter_pairs)):
            ats = np.array(ats); ls = np.array(ls); cs = np.array(cs)
            assert len(ats) == len(ls) == len(cs) == 2
            assert ls[0] in [1, 2] and ls[1] in [1,2]; assert 0 <= ats[0] < n_at and 0 <= ats[1] < n_at
            idxs = 3*ats + cs; y = np.array([f[idxs[0], idxs[1]] for f in self.force_matrices])
            lab = r'$%s^{(%d)}_%s \sim %s^{(%d)}_%s$'%(atomic_sites[ats[0]], ls[0], cl[0], atomic_sites[ats[1]], ls[1], cl[1])
            cf = fig.axes[i].tricontourf(bx, by, y, levels=levels, cmap="RdGy")
            if rows == 1:
                fig.colorbar(cf, ax=fig.axes[i], format='%.0e')
            fig.axes[i].set_aspect('equal')
            fig.axes[i].set_xlabel(r"$b_x$"); fig.axes[i].set_ylabel(r"$b_y$")
            fig.axes[i].text(0.7, 1.05, lab, transform=fig.axes[i].transAxes, size=10, weight='ultralight')
            fig.axes[i].text(0.05, 1.05, "%.0e to %.0e"%(min(y), max(y)), transform=fig.axes[i].transAxes, size=10, weight='ultralight')
        plt.tight_layout()
        fig.savefig(self.__out_dir + f"forces{'_' if pfx != '' else ''}{pfx}.png")
    
    # Do all available functions.
    def output_all_analysis(self, levels=DEFAULT_ABS_MIN_ENERGY):
        self.save_raw_data()
        self.plot_e_vs_b(levels=levels)
        self.plot_z_vs_b(levels=levels)
        self.plot_e_vs_z()


# Use basis linear regression to fit GSFE to leading 6 fourier terms
# See Eq(4), Carr 2018
# Performs a mapping to a BZ associated with 60-degree realspace lattice, which works regardless of 
# lattice in original POSCAR
# * Coefficients are in eV, while everything else is in meV above
class FourierGSFE:
    def __init__(self, energies, b_matrix, lattice_matrix, ltype='hexagonal', sampling_type='grid', dump=False):
        self.stype = sampling_type; assert isinstance(self.stype, str)
        self.__A = lattice_matrix
        print("Initializing FourierGSFE object")
        assert ltype == 'hexagonal', "Non-hexagonal lattices not supported"
        assert len(energies) == len(b_matrix), f"Must have same number of configurations as energies, but got {len(b_matrix)} vs. {len(energies)}"
        b_matrix = b_matrix[:,:2]
        if dump:
            print(f"Configurations (direct basis):\n {b_matrix}\nConfigurations (Cartesian):\n {(self.__A @ b_matrix.T).T}")
        minenergy = min(energies); self.GSFE = np.array(energies) - minenergy
        if dump:
            print(f"Adjusted energies (eV): {self.GSFE}")
        self.nb = len(b_matrix); self.b_matrix = b_matrix
        self.__fitted = False; self.coeffs = None; self.reg = None
        self.X = self.__b_to_fourier_basis(b_matrix)
    def __b_to_vw(self, b_mat):
        print(f"Scaling by lattice constant {LA.norm(self.__A[:,0])}")
        assert self.__A.shape == (2,2)
        M = 2 * pi * np.array([[1,-1/sqrt(3)],[0,2/sqrt(3)]]) / LA.norm(self.__A.T[0])
        return (M @ self.__A @ b_mat.T).T
    def __vw_to_fourier_basis(self, vw, sym=False):
        nbas = 3 if sym else 5
        X = np.zeros((len(vw), nbas)); v = vw[:,0]; w = vw[:,1]
        X[:,0] = np.cos(v) + np.cos(w) + np.cos(v + w)
        X[:,1] = np.cos(v + 2*w) + np.cos(v - w) + np.cos(2*v + w) 
        X[:,2] = np.cos(2*v) + np.cos(2*w) + np.cos(2*v + 2*w) 
        if not sym:
            X[:,3] = np.sin(v) + np.sin(w) - np.sin(v + w)
            X[:,4] = np.sin(2*v + 2*w) - np.sin(2*v) - np.sin(2*w)
        # X[:,5] = np.cos(3*v) + np.cos(3*w) + np.cos(3*v + 3*w)
        # X[:,6] = np.cos(v - 2*w) + np.cos(2*v + 3*w) + np.cos(3*v + w)
        # X[:,7] = np.cos(2*v - w) + np.cos(v + 3*w) + np.cos(3*v + 2*w)
        # X[:,8] = np.sin(3*v + 3*w) - np.sin(3*v) - np.sin(3*w)
        return X
    def __b_to_fourier_basis(self, b_mat):
        return self.__vw_to_fourier_basis(self.__b_to_vw(b_mat))
    def __sp_fit_coeff(self):
        assert not self.__fitted, f"Fitting already done"
        print("Fitting energies to Fourier series...")
        self.reg = LinearRegression().fit(self.X, self.GSFE)
        self.coeffs = np.append(self.reg.intercept_, self.reg.coef_)
        self.__fitted = True
        print(f"Final coefficients:\n{self.coeffs}")
        return self.coeffs
    def __biasify(self, X):
        npts = len(X); z = np.ones(npts).reshape((npts, 1))
        return np.concatenate((z, X), axis=1) # bias trick
    def __fit_coeff(self):
        self.__fitted = True
        X0 = self.__biasify(self.X)
        X = X0#[:,:-1]
        self.X = X0
        print(f"Fitting {len(X[0])} coefficients (1 less than total)")
        self.coeffs = np.dot(LA.pinv(np.dot(X.T, X)), np.dot(X.T, self.GSFE))
        #self.coeffs = np.append(self.coeffs, max(self.GSFE)/3 - sum(self.coeffs))
        print(f"Coefficients: {self.coeffs}")
        print(f"Residual sum of squares: {self.get_score()}")
    def __ensure_fitted(self):
        if not self.__fitted:
            self.__fit_coeff()
        assert self.coeffs is not None, f"Fitting failed"
    def save_raw_data(self, outdir):
        print("Saving coefficents and score to file...")
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
        self.__ensure_fitted()
        np.save(outdir + FGSFE_COEFF_NAME, self.coeffs)
        np.savetxt(outdir + FGSFE_COEFF_TXT, self.coeffs)
        np.save(outdir + FGSFE_SCORE_NAME, self.get_score())
    def get_coeffs(self):
        self.__ensure_fitted(); return self.coeffs
    def sp_get_score(self):
        self.__ensure_fitted()
        print(f"Fit score = {self.reg.score(self.X, self.GSFE)}")
        return self.reg.score(self.X, self.GSFE)
    def get_score(self):
        self.__ensure_fitted()
        return np.sum(np.square(self.GSFE - np.dot(self.X, self.coeffs)))
    def sp_predict(self, b_matrix):
        self.__ensure_fitted()
        return self.reg.predict(self.__b_to_fourier_basis(b_matrix[:,:2]))
    def predict(self, b_matrix):
        X = self.__biasify(self.__b_to_fourier_basis(b_matrix[:,:2]))
        return X @ self.coeffs
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
        print(f"Predicted GSFE:\n {x}\nActual GSFE: \n{y}")
        print(f"PRED (min={min(x)} max={max(x)} @ {np.where(x == max(x))[0]} range={max(x)-min(x)}) ACTUAL (min={min(y)} max={max(y)} @ {np.where(y == max(y))[0]} range={max(y)-min(y)})")
        fig.savefig(outdir + outname)
    def plot_percent_error(self, b_mat, actual, title='', outpath='percent_err.png', bins=10):
        pred = self.predict(b_mat)
        plt.clf(); plt.title(f"Percent error {'' if title == '' else 'in ' + title}")
        plt.xlabel('Percent error'); plt.ylabel('Frequency')
        diffs = np.array(pred) - np.array(actual)
        print(f"Pred-actual diffs:\n{diffs}")
        pred = pred[actual > 1e-2]; actual = actual[actual > 1e-2]
        plt.hist(100*abs((pred-actual)/(actual)), bins=bins)
        plt.savefig(outpath)
    def output_all_analysis(self, outdir):
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
        self.__ensure_fitted()
        self.save_raw_data(outdir)
        self.plot_pred_vs_actual(outdir)
        print("All Fourier GSFE analysis written to file")

