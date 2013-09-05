import numpy as np
import numpy.linalg as la
import find_neighbors as fn

class ChangingNeighbors:
    
    def __init__(self, file1, file2, cutoffs):
        # generate a list (size nr atoms) of arrays containing the
        # indicies of neighbors
        """
        Returns a numpy array with the first two cols defining the pair and the 
        three rest defining the cna number.
        
        @type  file1: path to a .con file
        @param file1: The FindNeighbors object for curresponding atoms system
        """
        self.neigh1 = fn.FindNeighbors(file1, cutoffs)
        self.neigh2 = fn.FindNeighbors(file2, cutoffs)
        self.sameNeighbors = self.getIntersects()
        self.lostNeighbors = self.getOldMinusNew()
        self.newNeighbors = self.getNewMinusOld()
        
    def getIntersects(self):
        res = []
        for i in xrange(self.neigh1.nrAtoms):
            #            res.append(np.intersect1d(self.neigh1.getAllNeighbors(i),
            #                          self.neigh2.getAllNeighbors(i)))
            res.append(np.intersect1d(self.neigh1.getIndiciesAllNeighborsForAtom(i),
                                      self.neigh2.getIndiciesAllNeighborsForAtom(i)))
        return res
    
    
    
    
    def getDiffs(self, neighSet1, neighSet2):
        res = []
        for i in xrange(neighSet1.nrAtoms):
            #    res.append(np.setdiff1d(neighSet1.getAllNeighbors(i),
            #                        neighSet2.getAllNeighbors(i)))
            res.append(np.setdiff1d(neighSet1.getIndiciesAllNeighborsForAtom(i),
                                    neighSet2.getIndiciesAllNeighborsForAtom(i)))
        return res
        
    def getOldMinusNew(self):
        return self.getDiffs(self.neigh1, self.neigh2)
        
    def getNewMinusOld(self):
        return self.getDiffs(self.neigh2, self.neigh1)
        
    # pass in one of the variables created in the constructor
    # to get the number of elements in the 1d array for each atom
    def nrValues(self, sets):
        return np.array(map(lambda x: x.size, sets))
