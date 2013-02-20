import os
import numpy as np
import numpy.linalg as la
import eon_tools as et
try:
    import converter as cc
except ImportError:
    import conf_convert as cc

class FindNeighbors:

    def __init__(self, atoms, cutoffs=None, distances=None, doSaveDistances=True):
        self.distances = np.zeros(0)
        self.fileName = None
        self.path = None
        if type(atoms) is et.Atoms:
            self.atoms = atoms
            self.nrAtoms = self.atoms.num_atoms
            if isinstance(distances, np.ndarray):
                self.distances = distances
            
        else:
            # we need to load from a file
            self.atoms = cc.con_2_eon_atoms(atoms)
            self.nrAtoms = self.atoms.num_atoms
            self.path, self.fileName = os.path.split(atoms)
            self.distName = self.fileName.split('.')[0] + '_distances.dat'
            if os.path.isfile(os.path.join(self.path, self.distName)):
                self.loadDistances(os.path.join(self.path, self.distName))
        self.cutoffs = cutoffs or {}.copy()
        
        self.atomTypes = np.unique(self.atoms.atomic_numbers)
        self.atomIndxs = self._getAtomIndexes()
        
        
        if self.distances.shape != (self.nrAtoms,self.nrAtoms):
            self.distances = np.zeros((self.nrAtoms, self.nrAtoms))
            self._calcDistances()
            if doSaveDistances and self.path:
                self.saveDistances(os.path.join(self.path, self.distName))
    
    def setCutOffs(self, cutoffs):
        """
        For the case H (1) and Aluminum (13).
        H-H 2.2
        H-Al 3.7
        Al-Al 4.1
            
        {1:  {1:2.2 , 13: 3.7},
         13: {1:3.7 , 13: 4.1}}
        """
        self.cutoffs = cutoffs
        
    def _getAtomIndexes(self):
        atomNrs = self.atoms.atomic_numbers
        indexArr = np.array(range(atomNrs.size))
        indxs = {}
        for i in np.unique(atomNrs):
            indxs[i] = indexArr[atomNrs == i]
        return indxs
        
    def _calcDistances(self):
        for i in xrange(self.nrAtoms):
            for j in xrange(i):
                self.distances[i][j] = \
                self.atoms.get_distance_atoms(i,j)
    
    def getNeighbors(self, i, nType):
        """
        Returns neighbors of nType from atoms i
        """
        aType = self.atoms.atomic_numbers[i]
        # Remember distances is a tringular matrix
        dists = np.hstack((self.distances[i,:i], self.distances[i:,i]))
        nrNeighbors = dists[dists < self.cutoffs[aType][nType]].size
        neighbors = np.argsort(dists)[:nrNeighbors]
        res = np.intersect1d(self.atomIndxs.get(nType,[]), neighbors)
        res = res[res != i]
        return res
        
    def getAllNeighbors(self, i):
        """
        Returns all neighbors for atom i
        """
        res = []
        for nType in self.atomTypes:
            res.append(self.getNeighbors(i, nType))
        return np.hstack(res)
        
        
    def countNeighborsByType(self, atomType, neighborTypes):
        """
        Determines the nr-distribution of neighbours with the given type.
        Output can be used to generate a histogram.
        """
        data = {}
        for i in self.atomIndxs[atomType]:
            neighb = np.zeros(0, dtype='int')
            for nType in neighborTypes:
                neighb = np.hstack((neighb, self.getNeighbors(i, nType)))
            nrNeighb = neighb.size
            data[nrNeighb] = data.get(nrNeighb,0) + 1
        dArray = np.zeros(max(data.keys())+1,dtype='int32')
        for key in data:
            dArray[key] = data[key]
        return dArray
        
    def rdf(self, rMax, nrBins, refAtomTypes, scannedAtomTypes):
        rAtoms = np.zeros(0, dtype='int')
        sAtoms = np.zeros(0, dtype='int')
        for i in refAtomTypes:
            rAtoms = np.hstack((rAtoms, self.atomIndxs[i]))
        for i in scannedAtomTypes:
            sAtoms = np.hstack((sAtoms, self.atomIndxs[i]))
        dr = float(rMax)/nrBins
        dr3 = dr**3
        
        bins = np.zeros(nrBins)
        for i in rAtoms:
            aBins = np.zeros(nrBins)
            dists = np.hstack((self.distances[i,:i], self.distances[i:,i]))[sAtoms]
            for j in xrange(nrBins):
                aBins[j] = dists[(dists > (dr*j)) * (dists < (dr*(j+1)))].size
            bins += aBins
        bins = bins / rAtoms.size
        for i in xrange(bins.size):
            bins[i] = bins[i] / (4/3 * np.pi * dr3 *((i+1)**3 - i**3))
                
        return bins
        
        
    def saveDistances(self, filePath):
        f = open(filePath, 'w')
        for i in xrange(self.nrAtoms):
            self.distances[i,:i].tofile(f)
        f.close()
        
    def loadDistances(self, filePath):
        self.distances = np.zeros((self.nrAtoms,self.nrAtoms))
        f = open(filePath, 'r')
        for i in xrange(self.nrAtoms):
            self.distances[i,:i] = np.fromfile(f, count=i)
        f.close()
        
