# Compute dynamical matrix block element for given q = Gtilde and phonopy object ph
def dm(q, ph):
    pass #! import

# Create level-1 block matrix
def block_l1():
    pass # TODO

# Create level-2 block matrix with intralayer and interlayer terms
def block_l2(D_intras, D_inter):
    assert len(D_intras) == 2
    assert D_intras[0].shape == D_intras[1].shape == D_inter.shape
    return np.block([[D_intras[0], D_inter], [D_inter.conjudate().T, D_intras[1]]])
