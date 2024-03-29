import numpy as np
import numpy.linalg as LA
import matplotlib.pyplot as plt
import os
import phonopy
from math import sqrt, pi
from sklearn.linear_model import LinearRegression
from ___constants_names import FFTENS_COEFF_NAME, FFTENS_COEFF_TXT, FFTENS_SCORE_NAME
from ___constants_output import NTERMS
from __directory_searchers import checkPath
from scipy.interpolate import griddata

# TODO rewrite using map() instead of brute force iteration

# Realspace interpolation with linear spline fitting
class ForceInterp:
    def __init__(self, ph_list, b_matrix):
        for i in range(len(ph_list)):
            ph_list[i].symmetrize_force_constants()
        self.force_matrices = [ph.get_dynamical_matrix_at_q([0,0,0]) for ph in ph_list]
        fcs = self.force_matrices[0].shape
        self.fc_tnsr = np.real_if_close([[[fc[j][k] for fc in self.force_matrices] for k in range(fcs[1])] for j in range(fcs[0])])
        print("Initializing spline force interpolator object")
        print(f"Force tensor shape: {self.fc_tnsr.shape}")
        self.b_matrix = b_matrix[:,:2]
        # self.f_mat = self.__functional_matrix()
    def __b_to_xy(self, b_matrix):
        x, y = b_matrix[:,0], b_matrix[:,1]
        return x, y
    # def __functional_matrix(self):
    #     x, y = self.__b_to_xy(self.b_matrix)
    #     return [[interp2d(x, y, fc) for fc in fcrow] for fcrow in self.fc_tnsr]
    def fc_tnsr_at(self, new_b_matrix):
        new_b_matrix = new_b_matrix[:,:2]
        def itp(frcs):
            interp = griddata(self.b_matrix, frcs, new_b_matrix, method='cubic')
            assert not np.any(np.isnan(interp)), f"NaN detected in interp {interp}"
        # fc_tnsr = np.array([[f(*self.__b_to_xy(b_matrix)) for f in frow] for frow in self.f_mat])
        fc_tnsr = np.array([[itp(frcs) for frcs in frow] for frow in self.fc_tnsr])
        print(f"Interpolated force tensor shape: {fc_tnsr.shape}")
        return fc_tnsr


# Fourier space interpolation with k-shell fitting
class FourierForceInterp:
    def __init__(self, ph_list, b_matrix, lattice_matrix, ltype='hexagonal', sampling_type='grid', dump=False):
        for i in range(len(ph_list)):
            ph_list[i].symmetrize_force_constants()
            ph_list[i].produce_force_constants()
        self.stype = sampling_type; assert isinstance(self.stype, str)
        self.__A = lattice_matrix
        print("Initializing Fourier force interpolator object")
        assert ltype == 'hexagonal', "Non-hexagonal lattices not supported"
        assert len(ph_list) == len(b_matrix), f"Must have same number of configurations as energies, but got {len(b_matrix)} vs. {len(ph_list)}"
        b_matrix = b_matrix[:,:2]
        if dump:
            print(f"Configurations (direct basis):\n {b_matrix}\nConfigurations (Cartesian):\n {(self.__A @ b_matrix.T).T}")
        
        # shape: (nS, nS, 3, 3, ncfg) with ncfg number of configs, nS num atoms in supercell
        self.fc_tnsr = np.transpose([ph.force_constants for ph in ph_list], axes=(1,2,3,4,0)) 
        self.shape = self.fc_tnsr.shape
        print(f"Force tensor shape: {self.shape}")
        self.nb = len(b_matrix); self.b_matrix = b_matrix
        self.__fitted = False; self.coeffs = None; self.reg_tnsr = None; self.score_mat = None
        self.X = self.__b_to_fourier_basis(b_matrix)
    def __b_to_vw(self, b_mat):
        print(f"Scaling by lattice constant {LA.norm(self.__A[:,0])}")
        assert self.__A.shape == (2,2)
        M = 2 * pi * np.array([[1,-1/sqrt(3)],[0,2/sqrt(3)]]) / LA.norm(self.__A.T[0])
        return (M @ self.__A @ b_mat.T).T
    def __vw_to_fourier_basis(self, vw, nsh=3):
        assert 1 <= nsh <= 3 and isinstance(nsh, int)
        v = vw[:,0]; w = vw[:,1]
        X = np.ones((len(vw), NTERMS[nsh]))
        X[:,0] = np.cos(v) + np.cos(w) + np.cos(v + w)
        X[:,1] = np.sin(v) + np.sin(w) - np.sin(v + w)
        if nsh >= 2:
            X[:,2] = np.cos(v + 2*w) + np.cos(v - w) + np.cos(2*v + w)
            X[:,3] = np.cos(2*v + 2*w) + np.cos(2*v) + np.cos(2*w)
            X[:,4] = np.sin(2*v + 2*w) - np.sin(2*v) - np.sin(2*w)
        if nsh >= 3:
            X[:,5] = np.cos(3*v) + np.cos(3*w) + np.cos(3*v + 3*w)
            X[:,6] = np.cos(v - 2*w) + np.cos(2*v + 3*w) + np.cos(3*v + w)
            X[:,7] = np.cos(2*v - w) + np.cos(v + 3*w) + np.cos(3*v + 2*w)
            X[:,8] = np.sin(3*v + 3*w) - np.sin(3*v) - np.sin(3*w)
        return X
    def __b_to_fourier_basis(self, b_mat):
        return self.__vw_to_fourier_basis(self.__b_to_vw(b_mat))
    def __fit_coeff(self):
        assert not self.__fitted, f"Fitting already done"
        print("Fitting energies to Fourier series...")
        self.reg_tnsr = [[[\
            [LinearRegression().fit(self.X, forces) for forces in f3] \
            for f3 in f2] \
                for f2 in f1] \
                    for f1 in self.fc_tnsr] # shape: (nS, nS, 3, 3)
        self.coeffs = np.array(\
            [[[[\
            np.append(reg.intercept_, reg.coef_) for reg in r3] \
                for r3 in r2] \
                    for r2 in r1] \
                        for r1 in self.reg_tnsr])
        print(f"Fourier force constants fit coefficient tensor shape: {self.coeffs.shape}")
        self.__fitted = True; self.score_mat = None
        return self.coeffs
    def __ensure_fitted(self):
        if not self.__fitted:
            self.__fit_coeff()
        assert self.coeffs is not None and self.reg_tnsr is not None, f"Fitting failed"
    def save_raw_data(self, outdir):
        print("Saving coefficents and score to file...")
        outdir = checkPath(outdir); assert os.path.isdir(outdir), f"Directory {outdir} does not exist"
        self.__ensure_fitted()
        np.save(outdir + FFTENS_COEFF_NAME, self.coeffs)
        np.savetxt(outdir + FFTENS_COEFF_TXT, self.coeffs)
        np.save(outdir + FFTENS_SCORE_NAME, self.get_score())
    def get_coeffs(self):
        self.__ensure_fitted(); return self.coeffs
    def get_score(self):
        self.__ensure_fitted()
        if self.score_mat is None:
            self.score_mat = np.array(\
                [[[[\
                reg.score(self.X, f) for f, reg in zip(f3, r3)] \
                    for f3, r3 in zip(f2, r2)] \
                        for f2, r2 in zip(f1, r1)] \
                            for f1, r1 in zip(self.fc_tnsr, self.reg_tnsr)])
            print("Constructed Fourier force fit score tensor")
        return self.score_mat
    def predict(self, b_matrix, dbglog=False):
        self.__ensure_fitted()
        X = self.__b_to_fourier_basis(b_matrix[:,:2])

        # Shape: (nS, nS, 3, 3, ncfg*) with ncfg* == len(b_matrix) in the function arg
        predtnsr = np.array(\
            [[[[\
            reg.predict(X) for reg in r3] \
                for r3 in r2] \
                    for r2 in r1] \
                        for r1 in self.reg_tnsr])
        if dbglog:
            print(f"Prediction tensor shape: {predtnsr.shape}")
        return predtnsr
