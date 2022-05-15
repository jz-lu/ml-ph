from pymatgen.io.vasp.outputs import Chgcar

c = Chgcar.from_file('CHGCAR')

c.write_file('TESTME')