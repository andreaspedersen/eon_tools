"""
Author: Kjartan Myrdal
Date: June 2012
License: Undefined, if you have this you're probably allowed to use it.
"""

import numpy as np
try:
    from ase.data import atomic_numbers
    import ase
    isAse = True
except ImportError:
    print 'Ase module not found, con_2_ase will not be defined.'
    isAse = False
import eon_tools as et

_floatErr = 1e-6


if isAse:
    """
    Loads con file
    Takes in a path to a .con file and returns
    ase.Atoms object. Optional pbc variable.
    """
    def con_2_ase(path, boundary = [1,1,1]):
        
        if not path.endswith('.con'):
            path += '.con'
        f = open(path, 'r')
        
        ## seed to random
        f.readline()

        ## time stuff
        f.readline()

        basis_vec = map(float, f.readline().split())

        angle_basis_vec = map(lambda x: x - 90 < _floatErr, map(float, f.readline().split()))
        if not np.array(angle_basis_vec).sum():
            print "Only orthogonal lattice vectores are accepted!"
            return

        f.readline()
        f.readline()
        nrComponents = int(f.readline().strip())
        nrAtomsPrType = np.array(map(int, f.readline().split()))

        #Masses, skiped in original parser. Same here
        f.readline()
        
        #All data is parced into these arrays before Atoms object is created.
        pos = np.zeros((nrAtomsPrType.sum(), 3))
        atomNumbers = np.zeros(nrAtomsPrType.sum())
        filt = np.zeros(nrAtomsPrType.sum(), dtype='int32')

        #
        for comp in xrange(nrComponents):
            element = f.readline().strip()
            if element.isdigit():
                element = int(element)
            else:
                element = atomic_numbers.get(element) or self._unkownElement(element)
            
            f.readline()
            
            for i in xrange(nrAtomsPrType[0:comp].sum(), nrAtomsPrType[0:comp+1].sum()):
                data = f.readline().split()
                pos[i] = map(float, data[:3])
                atomNumbers[i] = element
                filt[i] = int(data[3])
                
        f.close()
        return (ase.Atoms(numbers=atomNumbers, positions=pos, cell=basis_vec, pbc=boundary), abs(filt-1))

    # convert atoms object
    def eon_2_ase(e_atoms, boundary = [1,1,1]):

        return ase.Atoms(numbers=e_atoms.get_atomic_numbers(), positions=e_atoms.r.copy(), cell=e_atoms.box.copy(), pbc=boundary)

    # convert atoms object
    def ase_2_eon(ase_atoms, free = None):
        
        e_atoms = et.Atoms(0)
        e_atoms.r = ase_atoms.get_positions()
        e_atoms.num_atoms = ase_atoms.get_number_of_atoms()
        e_atoms.box = ase_atoms.get_cell()
        e_atoms.free = free or np.ones(e_atoms.num_atoms,dtype='int')
        e_atoms.atomic_numbers = ase_atoms.get_atomic_numbers()
        e_atoms.mass = ase_atoms.get_masses()
        
        return e_atoms
            


def _unkownElement(element):
    print "Element %s is unknown." %(element,)
    return 1
    
    
def con_2_eon_atoms(path):
    """
        Loads con file
        Takes in a path to a .con file and returns
        eon_atoms object.
    """
    at = et.Atoms(path)
#    at = et.Atoms(0)
#    at.read_con(path)
    return at
    

    
    
    

