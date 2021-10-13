import numpy as np
import os
from ___constants_names import (
    RELAX_CODE_PATH, RELAX_CODE_SUBPATH, 
    RELAX_CODE_OUT, RELAXED_DELTA_OUT, UNRELAXED_CONFIGS_OUT, 
    RELAXED_CONFIGS_NPY
)
from __directory_searchers import checkPath
from __class_Configuration import Configuration
import numpy.linalg as LA
import matplotlib.pyplot as plt
from pymatgen.io.vasp.inputs import Poscar
import argparse
from ___helpers_parsing import succ
from time import time

# RELAX_CODE_PATH = "/Users/jonathanlu/Documents/ml-ph/" + RELAX_CODE_SUBPATH

class RelaxerAPI:
    def __init__(self, theta, gridsz, outdir, cob):
        assert isinstance(gridsz, int) and gridsz > 1, f"Invalid grid size {gridsz}: must be positive integer"
        assert theta > 0 and theta < 180, f"Invalid twist angle {theta}"
        assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        self.pfx = f'_%.3lf'%theta
        self.outdir = checkPath(os.path.abspath(outdir))
        self.gridsz = gridsz; self.cob = cob; self.theta = theta
        self.langle = Configuration.lattice_angle_from_cob(self.cob)
        assert np.isclose(self.langle, 60) or np.isclose(self.langle, 120), f"Lattice angle must be 60 or 120, but is {self.langle}"
        print(f"Starting relaxer program {RELAX_CODE_PATH} in Julia...")
        stream = os.popen(f"julia {RELAX_CODE_PATH} -d {theta} -N {gridsz} -o {self.outdir}")
        print(stream.read())
        self.outpath = self.outdir + RELAX_CODE_OUT
        assert os.path.isfile(self.outpath), f"Failed to find expected relaxer output file at {self.outpath}"
        print("Relaxer code finished.")
        self.u = np.load(self.outdir + RELAXED_DELTA_OUT)
        self.b = np.load(self.outdir + UNRELAXED_CONFIGS_OUT)
        bprime_cart = np.load(self.outdir + RELAX_CODE_OUT)
        # Relaxer uses y-major order, we use x-major order, so convert
        bprime_cart = np.transpose(\
            bprime_cart.reshape((self.gridsz,self.gridsz,2)),\
                 axes=(1,0,2)).reshape((self.gridsz**2,2))
        self.bprime_cart = bprime_cart
        assert bprime_cart.shape == (self.gridsz**2, 2), f"Invalid relaxation matrix shape (expected {(self.gridsz**2, 2)}):\n {bprime_cart}"
        self.bprime_dir_raw = (LA.inv(self.cob) @ bprime_cart.T).T
        bprime = np.round((self.bprime_dir_raw + 1.0000001) % 1, 7) # mod unit cell torus
        self.bprime_dir = np.hstack((bprime, np.zeros(self.gridsz**2).reshape(self.gridsz**2, 1)))
        self.plot_quiver()
    def get_configs(self, save=True, cartesian=False):
        if save:
            np.save(self.outdir + RELAXED_CONFIGS_NPY, self.bprime_dir)
        return self.bprime_cart if cartesian else self.bprime_dir
    def plot_relaxation(self, filename='relax.png'):
        plt.clf(); _, ax = plt.subplots()
        plt.scatter(self.b[:,0], self.b[:,1], c='royalblue', alpha=0.15, label='before')
        plt.scatter(self.bprime_cart[:,0], self.bprime_cart[:,1], c='royalblue', label='after')
        ax.set_aspect('equal') # prevent stretching of space in plot
        pts = [[1/3, 1/3], [2/3, 2/3]]
        if np.isclose(self.langle, 120):
            pts = [[1/3, 2/3], [2/3, 1/3]]
        pts = np.array(pts); pts = (self.cob @ pts.T).T # make Cartesian
        labels = ['AB', 'BA']
        for (x, y), lab in zip(pts, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-5,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        plt.legend(); plt.title("Relaxer before-after for " + str(self.theta) + r"$^\circ$")
        ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$y$')
        filename = filename[:filename.index('.')] + self.pfx + filename[filename.index('.'):]
        outname = self.outdir + 'ba_' + filename
        plt.savefig(outname)
        succ("Successfully wrote relax Cartesian before-after plot out to " + outname)
    def plot_relaxation_direct(self, filename='relax_direct.png'):
        filename = filename[:filename.index('.')] + self.pfx + filename[filename.index('.'):]
        plt.clf(); _, ax = plt.subplots()
        b = (LA.inv(self.cob) @ self.b.T).T
        plt.scatter(b[:,0], b[:,1], c='royalblue', alpha=0.15, label='before')
        plt.scatter(self.bprime_dir[:,0], self.bprime_dir[:,1], c='royalblue', label='after')
        ax.set_aspect('equal') # prevent stretching of space in plot
        pts = [[1/3, 1/3], [2/3, 2/3]]
        if np.isclose(self.langle, 120):
            pts = [[1/3, 2/3], [2/3, 1/3]]
        labels = ['AB', 'BA']
        for (x, y), lab in zip(pts, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-5,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        plt.legend(); plt.title("Relaxer before-after for " + str(self.theta) + r"$^\circ$")
        ax.set_xlabel(r'$\mathbf{a}_1$'); ax.set_ylabel(r'$\mathbf{a}_2$')
        outname = self.outdir + filename
        plt.savefig(outname)
        succ("Successfully wrote relax Cartesian before-after plot out to " + outname)
    def plot_quiver(self, filename='quiver.png'):
        filename = filename[:filename.index('.')] + self.pfx + filename[filename.index('.'):]
        plt.clf(); _, ax = plt.subplots()
        plt.quiver(self.b[:,0], self.b[:,1], self.u[:,0], self.u[:,1])
        ax.set_aspect('equal') # prevent stretching of space in plot
        pts = [[1/3, 1/3], [2/3, 2/3]]
        if np.isclose(self.langle, 120):
            pts = [[1/3, 2/3], [2/3, 1/3]]
        pts = np.array(pts); pts = (self.cob @ pts.T).T # make Cartesian
        labels = ['AB', 'BA']
        for (x, y), lab in zip(pts, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-5,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        plt.title("Relaxation field for " + str(self.theta) + r"$^\circ$")
        ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$y$')
        outname = self.outdir + 'vecfld_' + filename
        plt.savefig(outname)
        succ("Successfully wrote relax Cartesian vector field plot out to " + outname)
    def plot_quiver_direct(self, filename='quiver_direct.png'):
        filename = filename[:filename.index('.')] + self.pfx + filename[filename.index('.'):]
        b = (LA.inv(self.cob) @ self.b.T).T; u = (LA.inv(self.cob) @ self.u.T).T
        plt.clf(); _, ax = plt.subplots()
        plt.quiver(b[:,0], b[:,1], u[:,0], u[:,1])
        ax.set_aspect('equal') # prevent stretching of space in plot
        pts = [[1/3, 1/3], [2/3, 2/3]]
        if np.isclose(self.langle, 120):
            pts = [[1/3, 2/3], [2/3, 1/3]]
        labels = ['AB', 'BA']
        for (x, y), lab in zip(pts, labels):
            plt.annotate(lab, (x, y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(-5,0), # distance from text to points (x,y)
                 ha='center') # horizontal alignment can be left, right or center
        plt.title("Relaxation field for " + str(self.theta) + r"$^\circ$")
        ax.set_xlabel(r'$\mathbf{a}_1$'); ax.set_ylabel(r'$\mathbf{a}_2$')
        outname = self.outdir + filename
        plt.savefig(outname)
        succ("Successfully wrote relax direct vector field plot out to " + outname)
    def plot_all(self):
        self.plot_relaxation(); self.plot_relaxation_direct()
        self.plot_quiver(); self.plot_quiver_direct()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Relax configurations in continuum model")
    parser.add_argument("angle", type=float, help="twist angle")
    parser.add_argument("gridsz", type=int, help="enter N for NxN sampling")
    parser.add_argument("-d", "--dir", type=str, help="output directory", default='.')
    parser.add_argument("-p", "--poscar", type=str, help="POSCAR file path", default="POSCAR")
    parser.add_argument("-c", "--cfg", help="get relaxed configs", action="store_true")
    args = parser.parse_args()
    assert 0 < args.angle < 180, f"Invalid twist angle {args.angle}"
    assert args.gridsz > 0, f"Invalid grid size {args.gridsz}"
    assert os.path.isdir(args.dir), f"Directory {args.dir} does not exist"
    assert os.path.isfile(args.poscar), f"File {args.poscar} does not exist"
    start_time = time()

    cob = Poscar.from_file(args.poscar).structure.lattice.matrix.T[:2,:2]
    r_api = RelaxerAPI(args.angle, args.gridsz, args.dir, cob)
    if args.cfg:
        r_api.get_configs()
    r_api.plot_all()
    succ("== Configuration Analyzer Complete (Took %.3lfs) =="%(time()-start_time))

    






