import numpy as np
import os
from ___constants_names import (
    RELAX_CODE_PATH, 
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

class RelaxerAPI:
    def __init__(self, theta, gridsz, outdir, cob):
        assert isinstance(gridsz, int) and gridsz > 1, f"Invalid grid size {gridsz}: must be positive integer"
        assert theta > 0 and theta < 180, f"Invalid twist angle {theta}"
        assert os.path.isdir(outdir), f"Invalid directory {outdir}"
        self.outdir = checkPath(os.path.abspath(outdir))
        self.gridsz = gridsz; self.cob = cob
        self.langle = Configuration.lattice_angle_from_cob(self.cob)
        assert np.isclose(self.langle, 60) or np.isclose(self.langle, 120), f"Lattice angle must be 60 or 120, but is {self.langle}"
        print("Starting relaxer program in Julia...")
        stream = os.popen(f"julia {RELAX_CODE_PATH} -d {theta} -N {gridsz} -o {self.outdir}")
        print(stream.read())
        self.outpath = self.outdir + RELAX_CODE_OUT
        assert os.path.isfile(self.outpath), f"Failed to find expected relaxer output file at {self.outpath}"
        print("Relaxer code finished.")
        self.u = np.load(self.outdir + RELAXED_DELTA_OUT)
        self.b = np.load(self.outdir + UNRELAXED_CONFIGS_OUT)
        bprime_cart = np.load(self.outdir + RELAX_CODE_OUT)
        assert bprime_cart.shape == (self.gridsz**2, 2), f"Invalid relaxation matrix shape (expected {(self.gridsz**2, 2)}):\n {bprime_cart}"
        bprime = np.round(((LA.inv(self.cob) @ bprime_cart.T).T + 1.0000001) % 1, 7) # mod unit cell torus
        self.bprime = np.hstack((bprime, np.zeros(self.gridsz**2).reshape(self.gridsz**2, 1)))
    def get_configs(self, save=True):
        if save:
            np.save(self.outdir + RELAXED_CONFIGS_NPY, self.bprime)
        return self.bprime
    def plot_relaxation(self, filename='relax.png'):
        plt.clf(); _, ax = plt.subplots()
        plt.scatter(self.b[:,0], self.b[:,1], c='black', alpha=0.25, label='before')
        plt.scatter(self.bprime[:,0], self.bprime[:,1], c='royalblue', label='after')
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
        plt.legend(); plt.title("Relaxer before-after")
        ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$y$')
        outname = self.outdir + 'ba_' + filename
        plt.savefig(outname)
        succ("Successfully wrote relax before-after plot out to " + outname)
    def plot_quiver(self, filename='quiver.png'):
        plt.clf(); _, ax = plt.subplots()
        plt.quiver(self.b, self.u[:,0], self.u[:,1])
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
        plt.legend(); plt.title("Relaxer vector field")
        ax.set_xlabel(r'$x$'); ax.set_ylabel(r'$y$')
        outname = self.outdir + 'vecfld_' + filename
        plt.savefig(outname)
        succ("Successfully wrote relax before-after plot out to " + outname)
    def plot_all(self):
        self.plot_relaxation(); self.plot_quiver()


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

    






