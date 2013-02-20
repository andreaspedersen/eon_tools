import numpy as np

colors = {'gray': 1, 'green': 4, 'dark gray': 6, 'blue':7, 'red':8, 'yellow': 9}

def color(atoms,indexes=None,color='red', makeCopy=True):
    if makeCopy:
        atoms = atoms.copy()
    if indexes == None:
        indexes = np.array(range(atoms.get_number_atoms()))
    atoms.numbers[indexes] = colors[color]
    return atoms
