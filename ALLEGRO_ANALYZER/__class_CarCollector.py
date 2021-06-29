from ____exit_with_error import exit_with_error
from ___constants_names import POSCAR_UNIT_RELAXATION_NAME, POTCAR_NAME, KPOINTS_LINE_NAME, KPOINTS_MESH_NAME, INCAR_RELAXATION_NAME, POSCAR_NAME, CONTCAR_NAME
from ___constants_vasp import (
    RELAXATION_GRID_DENSITY, 
    RELAXATION_GRID_SHIFT, 
    NONRELAXATION_GRID_DENSITY, 
    NONRELAXATION_GRID_SHIFT, 
    INCAR_RELAX_SETTINGS, 
    INCAR_NORELAX_SCON_SETTINGS, 
    INCAR_VDW_SETTINGS, 
    POT_PMG_INIT_CMD
)
from ___constants_misc import ERR_NO_POSCAR
from __directory_searchers import checkPath
from __query_inputs import getInputName
from pymatgen.io.vasp.inputs import Poscar, Kpoints, Potcar, Incar, VaspInput
import os

# Not a typo, just a naming convention for "*CAR" of VASP I/O files.
class CarCollector:
    def __init__(self, ROOT, do_vdW, kpoints_is_gamma, need_line_kpoints, poscar=None):
        self.ROOT = checkPath(ROOT)
        self.do_vdW = do_vdW
        self.need_line_kpoints = need_line_kpoints
        self.kpoints_is_gamma = kpoints_is_gamma

        if poscar == None:
            print('Importing POSCAR from specified file directory %s...'%(ROOT))
            try:
                self.poscar = Poscar.from_file(self.ROOT + POSCAR_UNIT_RELAXATION_NAME)
            except Exception as err:
                print('Error:', err)
                exit_with_error(ERR_NO_POSCAR)
        else:
            self.poscar = poscar
        self.poscar.structure.sort() # sort poscar

        # Get the name of the crystal for naming purposes
        self.mat_name = getInputName(self.poscar)

    # Returns pymatgen Incar object
    def get_relaxation_incar(self, writeOut=False):
        # First let's build an array of default keys, which is helpful later
        default_keys = list(INCAR_RELAX_SETTINGS.keys())
        default_values = list(INCAR_RELAX_SETTINGS.values())
        vdW_keys = list(INCAR_VDW_SETTINGS.keys())
        vdW_values = list(INCAR_VDW_SETTINGS.values())

        incarExists = os.path.isfile(self.ROOT + INCAR_RELAXATION_NAME) # Look for existing INCAR
        if not incarExists:
            # Create a default
            print('No INCAR input found. Generating default relaxation settings...')
            incar = Incar.from_string('')
            for i in range(0, len(default_keys)):
                incar[default_keys[i]] = default_values[i]
        else:
            print('INCAR input found. Using parameters given and adding any other necessary ones...')
            try:
                incar = Incar.from_file(self.ROOT + INCAR_RELAXATION_NAME)
            except Exception as error:
                print('Error:', error)
                exit_with_error('Suggestion: you probably have a file labeled INCAR in the specified directory that is invalid.')

            for i in range(0, len(default_keys)):
                if default_keys[i] not in incar:
                    print('{} not in INCAR. Adding it automatically for relaxation calculations...'.format(default_keys[i]))
                    incar[default_keys[i]] = default_values[i]
        
        if self.do_vdW:
            print('Adding van der Waals interaction parameters to INCAR...')
            for i in range(0, len(vdW_keys)):
                incar[vdW_keys[i]] = vdW_values[i]
        
        # Get a SYSTEM tag for the name, as a comment
        incar['SYSTEM'] = self.mat_name

        if writeOut:
            print('Writing INCAR to file...')
            incar.write_file(self.ROOT + INCAR_RELAXATION_NAME)

        self.relaxation_incar = incar        
        return incar

    def get_norelax_scon_incar(self, writeOut=False):
        default_keys = list(INCAR_NORELAX_SCON_SETTINGS.keys())
        default_values = list(INCAR_NORELAX_SCON_SETTINGS.values())
        vdW_keys = list(INCAR_VDW_SETTINGS.keys())
        vdW_values = list(INCAR_VDW_SETTINGS.values())

        incarExists = os.path.isfile(self.ROOT + INCAR_RELAXATION_NAME) # Look for existing INCAR
        if not incarExists:
            # Create a default
            print('No INCAR input found. Generating default relaxation settings...')
            incar = Incar.from_string('')
            for i in range(0, len(default_keys)):
                incar[default_keys[i]] = default_values[i]
        else:
            print('INCAR input found. Using parameters given and adding any other necessary ones...')
            try:
                incar = Incar.from_file(self.ROOT + INCAR_RELAXATION_NAME)
            except Exception as error:
                print('Error:', error)
                exit_with_error('Suggestion: you probably have a file labeled INCAR in the specified directory that is invalid.')

            for i in range(0, len(default_keys)):
                if default_keys[i] not in incar:
                    print('{} not in INCAR. Adding it automatically for relaxation calculations...'.format(default_keys[i]))
                    incar[default_keys[i]] = default_values[i]
        
        if self.do_vdW:
            print('Adding van der Waals interaction parameters to INCAR...')
            for i in range(0, len(vdW_keys)):
                incar[vdW_keys[i]] = vdW_values[i]
        
        # Get a SYSTEM tag for the name, as a comment
        incar['SYSTEM'] = self.mat_name

        if writeOut:
            print('Writing INCAR to file...')
            incar.write_file(self.ROOT + INCAR_RELAXATION_NAME)

        self.relaxation_incar = incar        
        return incar

    # Returns pymatgen Kpoints object, by default we use Monkhorst-Pack
    def get_mesh_kpoints(self, grid=RELAXATION_GRID_DENSITY, shift=RELAXATION_GRID_SHIFT, writeOut=False):
        kpointsExists = os.path.isfile(self.ROOT + KPOINTS_MESH_NAME) # Look for existing KPOINTS
        if kpointsExists:
            try:
                kpoints = Kpoints.from_file(self.ROOT + KPOINTS_MESH_NAME)
                print('KPOINTS file inputted by user. Importing input for computation...')
            except Exception as err:
                print('Error:', err)
                exit_with_error('Suggestion: your user-inputted KPOINTS file is likely invalid. Check it.')

        else:
            if grid == RELAXATION_GRID_DENSITY and shift == RELAXATION_GRID_SHIFT:
                print('No KPOINTS file found. Generating default relaxation mesh...')
            else:
                print(f'No KPOINTS file found. Generating mesh according to user inputted grid {grid} and shift {shift}...')

            if self.kpoints_is_gamma:
                print('Using Gamma-centered scheme...')
                kpoints = Kpoints.gamma_automatic(grid, shift)
            else:
                print('Using Monkhorst-Pack scheme...')
                kpoints = Kpoints.monkhorst_automatic(grid, shift)
        
        kpoints.comment = self.mat_name # Standardize labeling style

        if (not kpointsExists) and writeOut:
            print('Writing new KPOINTS file...')
            kpoints.write_file(self.ROOT + KPOINTS_MESH_NAME)
        
        if not writeOut:
            print('KPOINTS object created without file writing.')

        self.kpoints_mesh = kpoints
        return kpoints
    
    # Returns pymatgen Potcar object fetched from Kaxiras lab group shared base
    def get_potcar(self, writeOut=False):
        if os.path.isfile(self.ROOT + POTCAR_NAME):
            potcar = Potcar.from_file(self.ROOT + POTCAR_NAME)
        else:
            try:
                atoms = self.poscar.site_symbols # Get array of atoms for potential, in order as given in POSCAR
            except Exception as err:
                print('Error:', err)
                exit_with_error('Possible source of error: your POSCAR is likely wrong. You cannot import POTCAR until POSCAR is made properly for the PUC. Check POSCAR input and POTCAR maker function input.')
            #print(atoms)

            # Get a pymatgen Potcar object with our lab group settings in .pmgrc.yaml
            func = 'PBE' if self.do_vdW else 'LDA_54'
            print('POTCAR SETTINGS: using ' + func + ' functional')
            try:
                potcar = Potcar(atoms, functional=func)
                if writeOut:
                    print('POTCARs fetched by pymatgen. Writing to file path %s'%(self.ROOT + POTCAR_NAME))
                    potcar.write_file(self.ROOT + POTCAR_NAME)
                else:
                    print('POTCARs fetched by pymatgen. Returned without writing to file.')
            except ValueError as err:
                print('Error:', err)
                exit_with_error('\nPossible solution: run the following command as a one-time intialization to set up directories, if using Odyssey. \n\n\t{}\n'.format(POT_PMG_INIT_CMD))

        self.potcar = potcar
        return potcar

    # Class method of building a full object for relaxation
    def build_relaxation_input(self, print_key_info=False, print_all_info=False):
        incar = self.get_relaxation_incar()
        kpoints = self.get_mesh_kpoints()
        potcar = self.get_potcar()
        v = VaspInput(incar, kpoints, self.poscar, potcar)

        if print_all_info:
            print('Relaxation Incar:', incar)
            print('Relaxation Kpoints mesh:', kpoints)
            print('Relaxation POSCAR, initial:', self.poscar)
        elif print_key_info:
            print('Relaxation POSCAR, intiial:', self.poscar)

        return v

    def build_norelax_input(self, print_key_info=False, print_all_info=False, 
                            grid=NONRELAXATION_GRID_DENSITY, shift=NONRELAXATION_GRID_SHIFT):
        incar = self.get_norelax_scon_incar()
        kpoints = self.get_mesh_kpoints(grid=grid, shift=shift)
        potcar = self.get_potcar()
        v = VaspInput(incar, kpoints, self.poscar, potcar)
        if print_all_info:
            print('Relaxation Incar:', incar)
            print('Relaxation Kpoints mesh:', kpoints)
            print('Relaxation POSCAR, initial:', self.poscar)
        elif print_key_info:
            print('Relaxation POSCAR, intiial:', self.poscar)
        return v

    # A static method to build a vasp input object
    @staticmethod
    def build_general_vasp_input_instance(incar, kpoints, poscar, potcar):
        return VaspInput(incar, kpoints, poscar, potcar)

    # Get the line kpoints if necessary, as determined by InputData class instance
    @staticmethod
    def get_line_kpoints(ROOT, need_line_kpoints):
        if need_line_kpoints:
            try:
                print('Parsing imported KPOINTS line file for band structure calculations...')
                kpoints_line = Kpoints.from_file(ROOT + KPOINTS_LINE_NAME)
            except Exception as err:
                exit_with_error('Error in importing line KPOINTS file for requested electronic and/or phononic band structure calculations: ' + err)
            return kpoints_line
        else:
            return None

    # Get a poscar object from file
    @staticmethod
    def poscar_from_file(dirName):
        return Poscar.from_file(checkPath(dirName) + POSCAR_NAME)

    # Get the layer distance
    @staticmethod
    def get_interlayer_spacing(DIR_RELAXATION):
        p = Poscar.from_file(checkPath(DIR_RELAXATION) + CONTCAR_NAME)
        atoms_z_pos = []
        atoms = p.structure.frac_coords

        for i in atoms:
            if i[2] > 0:
                atoms_z_pos.append(round(i[2], 6))
        atoms_z_pos = list(dict.fromkeys(atoms_z_pos))

        interlayer_spacing = min(atoms_z_pos)
        return interlayer_spacing