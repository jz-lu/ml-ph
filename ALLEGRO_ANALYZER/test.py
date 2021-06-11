# * Exporting/importing numpy data files 
import numpy as np

THIS_DIR = './testfeatures/'

# Some example data
a = []
a.append((np.asarray([0., 0., 0.]), 0.226122, -36.9177224))
a.append((np.asarray([0.      , 0.333333, 0.      ]), 0.22206, -36.9179183))
a.append((np.asarray([0.      , 0.666667, 0.      ]), 0.221466, -36.91792992))
a.append((np.asarray([0.333333, 0.      , 0.      ]), 0.221393, -36.91795228))
a.append((np.asarray([0.333333, 0.333333, 0.      ]), 0.21728, -36.9179966))
a.append((np.asarray([0.333333, 0.666667, 0.      ]), 0.223795, -36.91788055))
a.append((np.asarray([0.666667, 0.      , 0.      ]), 0.223927, -36.91789084))
a.append((np.asarray([0.666667, 0.333333, 0.      ]), 0.223875, -36.91787807))
a.append((np.asarray([0.666667, 0.666667, 0.      ]), 0.21511, -36.91800157))
THIS_FILE = THIS_DIR + 'bze-test'
np.save(THIS_FILE, np.array(a))
print('SAVED np-array to ' + THIS_FILE)
print(a)

# Try to restore it
b = np.load(THIS_FILE + '.npy', allow_pickle=True)
print('RECOVERED np-array:')
print(b)
assert np.all(a == b)
