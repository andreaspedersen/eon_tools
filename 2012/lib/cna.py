import numpy as np

def full_cna(fn):
    """
    Returns a numpy array with the first two cols defining the pair and the 
    three rest defining the cna number.
    
    @type  fn: FindNeighbors
    @param fn: The FindNeighbors object for curresponding atoms system
    """
    safn = []
    for i in xrange(fn.nrAtoms):
        neigh = fn.getAllNeighbors(i)
        neigh = neigh[neigh > i]
        cnas = np.empty((neigh.size, 5), dtype='int')
        cnas[:,0] = i
        cnas[:,1] = neigh
        for j in xrange(neigh.size):
            cnas[j,2:] = cna(i,neigh[j],fn)
        safn.append(cnas)
    return np.vstack(safn)

def cna(i, j, fn):
    """
    Calculates CNA for an individual pair of index i and j.
    
    @type   i: int
    @param  i: index of first atom
    @type   j: int
    @param  j: index of second atom
    @type  fn: FindNeighbors
    @param fn: The FindNeighbors object for curresponding atoms system
    """
    cNeigh = np.intersect1d(fn.getAllNeighbors(i),fn.getAllNeighbors(j))
    if not cNeigh.size:
        return 0,0,0
    bonds = {}
    cHeads = [cNeigh[0]]
    l = 1
    for k in xrange(cNeigh.size):
        if l == k:
            l += 1
            cHeads.append(cNeigh[k])
        bonds[cNeigh[k]] = np.intersect1d(fn.getAllNeighbors(cNeigh[k]), cNeigh[k+1:])
        queue = np.setdiff1d(bonds[cNeigh[k]],cNeigh[:l])
        tail = np.setdiff1d(cNeigh[l:],queue)
        cNeigh[l:] = np.hstack((queue,tail))
        l += queue.size
    
    bLengths = np.empty(len(cHeads))        
    for k in xrange(len(cHeads)):
        bLengths[k] = chainLength(cHeads[k], bonds).sum()
    
    nrBonds = 0
    for bond in bonds.values():
        nrBonds += bond.size
        
    return cNeigh.size, nrBonds, bLengths.max()
    
    
def chainLength(i, bonds):
    """
    Helperfunction for cna() to calculate the longest connected chain.
    Takes in index and a dictionary containing the connections and returns
    the two longest chians out from the given indexed atom.
    
    @type      i: int
    @param     i: index of the atom to scan out from
    @type  bonds: {int:numpy.array(dtype='int')}
    @param bonds: sort of a tree order for bonds amongst the atoms
    """
    b = bonds[i]
    if not b.size:
        return np.zeros(2)
    else:
        c = np.empty(b.size)
        for j in xrange(b.size):
            c[j] = chainLength(b[j],bonds)[0] + 1
        return np.sort(c)[:-3:-1]
        
def get_cna_by_index(i, cna):
    """
    Input array returned by full_cna
    Returns the CNA values for the base pairs to which i belongs    
    """
    return cna[(cna[:,0] == i) + (cna[:,1] == i)]
    
def get_cna_by_indexes(indexes, cna):
    vals = np.zeros(len(cna), dtype='bool')
    for i in indexes:
        vals += (cna[:,0] == i) + (cna[:,1] == i)
    return cna[vals]
    
def get_unique_cna(cna):
    """
    Input array returned by full_cna
    Returns counts for each occurring CNA value
        
    """
    un = dict()
    for row in cna[:,2:]:
        row = tuple(row)
        if row in un:
            un[row] += 1
        else:
            un[row] = 1
    for key in un:
        un[key] = float(un[key])/len(cna)
    return un
    
def cna_value_to_pair(value, cna):
    """
    Input array returned by full_cna
    Returns all base pairs with given value
    """
    return cna[:,:2][np.apply_along_axis(lambda x: x.all(),1,cna[:,2:] == value)]

def get_cna_count_for_value(value, cna):
    """
    Input array returned by full_cna
    Returns count for all base pairs with given value
    """
    
    count = []
    i = 0
    while 1:
        cna_i = get_cna_by_index(i, cna)
        if len(cna_i):
            count.append(len(cna_value_to_pair(value, cna_i)))
            i = i + 1
        else:
            break
    return np.array(count)

def get_index_fcc(cna):
    fcc = get_cna_count_for_value((4,2,1), cna)
    return [fcc[0] for fcc in enumerate(fcc) if fcc[1]==12]

def get_index_ico(cna):
    ico = get_cna_count_for_value((5,5,5), cna)
    return [ico[0] for ico in enumerate(ico) if ico[1]==12]

def get_index_almost_ico(cna, treshold = 12):
    ico = get_cna_count_for_value((5,5,5), cna) + get_cna_count_for_value((5,4,4), cna)
    return [ico[0] for ico in enumerate(ico) if ico[1]>=treshold]

def get_index_hcp(cna):
    cna_422 = get_cna_count_for_value((4,2,2), cna)
    cna_421 = get_cna_count_for_value((4,2,1), cna)
    
    index_422 = [cna_422[0] for cna_422 in enumerate(cna_422) if cna_422[1]==6]
    index_421 = [cna_421[0] for cna_421 in enumerate(cna_421) if cna_421[1]==6]
    
    return np.intersect1d(index_422, index_421)

def get_index_less_than_neighbors(cna, treshold = 8):
    indicies = []
    for i in range(len(cna)):
        neighbors = len(get_cna_by_index(i, cna))
        if neighbors > 0 and neighbors < treshold:
            indicies.append(i)
    return indicies

def cna_intersect(cna1, cna2):
    nrows, ncols = cna1.shape
    
    v = {'names':['f{0}'.format(i) for i in range(ncols)],
         'formats':ncols * [cna1.dtype]}
    
    res = np.intersect1d(cna1.view(v), cna2.view(v))
    
    return res.view(cna1.dtype).reshape(-1, ncols)
    
def cna_atoms_no_change(cna1, cna2):
    b1 = np.ones(len(cna1),dtype='bool')
    b2 = np.ones(len(cna2),dtype='bool')
    m = cna1[:,:1].max()
    
    # define new data structure (the five elements stored per based pair)
    v = {'names':['f{0}'.format(i) for i in range(5)],
         'formats':5 * [cna1.dtype]}
    # 'transform' old data to new data structure
    cna1v = cna1.view(v)
    cna2v = cna2.view(v)

    for i in xrange(len(cna1v)):
        b1[i] = cna1v[i] not in cna2v
        
    for i in xrange(len(cna2v)):
        b2[i] = cna2v[i] not in cna1v
    
    # the changed cna values
    changed = np.vstack((cna1[b1],cna2[b2]))
    
    # Indexes of atoms with changed cna values (not unique)
    indexes = np.hstack((changed[:,0],changed[:,1]))
    
    # Indexes of atoms with all same cna numbers (all minus changed)
    res = np.setdiff1d(np.arange(m),indexes)
    
    return res


