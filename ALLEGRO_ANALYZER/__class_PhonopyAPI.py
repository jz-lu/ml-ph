from ___constants_names import (
    SPOSCAR_NAME, PH_FORCE_SETS_NAME, 
    ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME, 
    MONOLAYER_DIR_NAME, CONFIG_DIR_NAME, CONFIG_SUBDIR_NAME
)
from ___constants_misc import ERR_PH_FORCE_SETS_NOT_FOUND
from ___helpers_parsing import greet, succ, warn, err
from __directory_searchers import checkPath, findDirsinDir
from __dirModifications import build_dir
import os, phonopy

class PhonopyAPI:
    def __init__(self, ROOT, spname=SPOSCAR_NAME):
        greet("Initialized Phonopy API object. Loading phonopy objects...")
        self.SPOSCAR_ERR_MSG = f"No supercell POSCAR with name {spname} found"
        self.ROOT = ROOT; self.spname = SPOSCAR_NAME
        self.nlayer, self.intra_list = self.__load_intra_ph_list(ROOT, spname)
        self.nconfig, self.inter_list = self.__load_inter_ph_list(ROOT, spname)
        succ("Successfully loaded interlayer and intralayer objects")

    # Load monolayer phonopy objects
    def __load_intra_ph_list(self, ROOT, spname=SPOSCAR_NAME):
        ROOT = checkPath(ROOT); ph_list = []
        layers = sorted(findDirsinDir(ROOT, MONOLAYER_DIR_NAME, searchType='start'))
        assert len(layers) > 1, "Twist calculations require at least 2 layers"
        assert len(layers) <= 2, "Twist calculations for more than 2 layers not supported (yet)"
        for layer in layers:
            ph_dir = build_dir([ROOT, layer, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
            assert os.path.isdir(ph_dir), f"Directory {ph_dir} not found"
            print(f"Extracting {spname} and {PH_FORCE_SETS_NAME} from {ph_dir}...")
            os.chdir(ph_dir)
            assert os.path.isfile(spname), self.SPOSCAR_ERR_MSG
            assert os.path.isfile(PH_FORCE_SETS_NAME), ERR_PH_FORCE_SETS_NOT_FOUND
            ph_list.append(phonopy.load(supercell_filename=spname))
        return len(layers), ph_list

    # Load configuration phonopy objects
    def __load_inter_ph_list(self, ROOT, spname=SPOSCAR_NAME):
        ROOT = build_dir([ROOT, CONFIG_DIR_NAME]); ph_list = []
        configs = sorted(findDirsinDir(ROOT, CONFIG_SUBDIR_NAME, searchType='start'))
        for config in configs:
            ph_dir = build_dir([ROOT, config, ANALYSIS_DIR_NAME, PHONOPY_DIR_NAME])
            assert os.path.isdir(ph_dir), f"Directory {ph_dir} not found"
            print(f"Extracting {spname} and {PH_FORCE_SETS_NAME} from {ph_dir}...")
            os.chdir(ph_dir)
            assert os.path.isfile(spname), self.SPOSCAR_ERR_MSG
            assert os.path.isfile(PH_FORCE_SETS_NAME), ERR_PH_FORCE_SETS_NOT_FOUND
            ph_list.append(phonopy.load(supercell_filename=spname))
        return len(configs), ph_list
    
    def nlayers(self):
        return self.nlayer
    def nconfigs(self):
        return self.nconfig
    def intra_ph_list(self):
        return self.intra_list
    def inter_ph_list(self):
        return self.inter_list

