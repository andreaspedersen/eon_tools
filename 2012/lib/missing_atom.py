import numpy as np
import numpy.linalg as la
try:
    import converter as cc
except:
    import conf_convert as cc
from ase import Atoms, Atom, view

class MissingAtom:

    def __init__(self, atoms, bound=3.3, nrNeighbors=12, atomTypes = [13]):
        if type(atoms) is Atoms:
            self.atoms = atoms
        else:
            self.atoms = cc.con_2_ase(atoms)[0]
        self.bound = bound
        self.nrNeighbors = nrNeighbors
        self.atomTypes = atomTypes
        self.atomIndxs = self._getAtomIndexes()
        self.nrAtoms = len(self.atomIndxs)
        self.distances = np.zeros((self.nrAtoms, self.nrAtoms))
        self.neighbors = np.zeros(self.nrAtoms, dtype='int32')
        self._calcDistances()
        self._countNeighbors(bound)
        self.missingNeighbors = self.indexOfNeighbMissing(nrNeighbors)
    
    def _getAtomIndexes(self):
        atomNrs = self.atoms.get_atomic_numbers()
        indxs = []
        for i in xrange(atomNrs.size):
            if atomNrs[i] in self.atomTypes:
                indxs.append(i)
        return np.array(indxs)
        
    def _calcDistances(self):
        for i in xrange(self.nrAtoms):
            for j in xrange(i):
                self.distances[i][j] = \
                self.atoms.get_distance(self.atomIndxs[i],
                                        self.atomIndxs[j], True)
                
    def _countNeighbors(self, bound):
        for i in xrange(self.nrAtoms):
            dists = np.hstack((self.distances[i,:i], self.distances[i+1:,i]))
            self.neighbors[i] = dists[dists < bound].size
            
    def indexOfNeighbors(self, i, bound):
        dists = np.hstack((self.distances[i,:i], self.distances[i+1:,i]))
        return np.argsort(dists)[:self.neighbors[i]]
        
    def getNeighborsOfPoint(self, point, bound):
        locations = self.atoms.get_positions()
        dists = np.zeros(self.nrAtoms)
        for i in xrange(self.nrAtoms):
            dists[i] = la.norm(self.particleToParticleVec(point,
                               locations[self.atomIndxs[i]]))
        return self.atomIndxs[dists < bound]
            
    def indexOfNeighbMissing(self, nrNeighbors):
        return self.atomIndxs[self.neighbors < nrNeighbors]
    
    def initialGuess(self, i):
        locations = self.atoms.get_positions()
        sol = np.zeros(3)
        neighbors = self.indexOfNeighbors(i, self.bound)
        vSum = np.zeros(3)
        for n in neighbors:
            vSum += self.particleToParticleVec(locations[i],locations[n])
        sol -= vSum/la.norm(vSum) * \
               np.hstack((self.distances[i,:i],
                          self.distances[i+1:,i])).min() - locations[i]
        return sol 
        
    def calcMeanLocation(self):
        bestGuess = None
        bestQuality = 0
        for i in self.missingNeighbors:
            guess = self.initialGuess(i)
            quality = self.guessQuality(guess)
            if quality > bestQuality:
                bestGuess = guess
                bestQuality = quality
        locations = self.atoms.get_positions()[self.missingNeighbors]
        for i in xrange(len(locations)):
            locations[i] = guess + \
                           self.particleToParticleVec(guess, locations[i])
        nrAtoms = len(self.missingNeighbors)
        x = locations[:,0].sum() / nrAtoms
        y = locations[:,1].sum() / nrAtoms
        z = locations[:,2].sum() / nrAtoms
        return np.array([x,y,z])
    
    def isGuessGood(self, guess):
        neigh = self.getNeighborsOfPoint(guess, self.bound)
        if neigh == self.missingNeighbors:
            return True
        return False
        
    def guessQuality(self, guess):
        neigh = self.getNeighborsOfPoint(guess, self.bound)
        return np.intersect1d(neigh, self.missingNeighbors).size

    
    def calcMeanLocation_new(self, atomsTemp=None, atomsToPop=1):
        trans_vac = np.zeros(3,'f')

        if atomsTemp == None:
            atomsTemp = Atoms(self.atoms)
        else:
            atomsTemp = Atoms(atomsTemp)
                            
        atomsTempOrg = Atoms(atomsTemp)
                
        for i in range(atomsToPop):
            atomsTemp.pop(-1)
            
#        box_x = self.atoms.get_cell()[0][0]
#        box_y = self.atoms.get_cell()[1][1]
#        box_z = self.atoms.get_cell()[2][2]
        pos_vac = None
        step = 0.001
        move_x = 0
        move_y = 0
        move_z = 0
        
        tries = 10000
        while tries:
            tries = tries - 1 
            ok = 1
            atomsTemp.translate([step * move_x, step * move_y, step * move_z])
            atomsTemp.set_scaled_positions(atomsTemp.get_scaled_positions())
            posTrans = atomsTemp.get_positions()#[self.missingNeighbors]
            
            #print len(atomsTemp)
            #print trans_vac
            x_min = min(posTrans[:,0])
            y_min = min(posTrans[:,1])
            z_min = min(posTrans[:,2])
            
            #print trans_vac
            #print [x_min,y_min,z_min]
                
            #if x_min < 1.5: #box_x / 3.:
            if (len(self.atomsInLayer(atomsTemp, 0, 0., 0.25)) != 8) or (len(self.atomsInLayer(atomsTemp, 0, 6.25, 0.25)) != 8):
                trans_vac[0] = trans_vac[0] + step
                ok = 0
                move_x = 1
            else:
                move_x = 0
        
            #if y_min < 1.5: #box_y / 3.:
            if (len(self.atomsInLayer(atomsTemp, 1, 0., 0.25)) != 8) or (len(self.atomsInLayer(atomsTemp, 1, 6.25, 0.25)) != 8):
                trans_vac[1] = trans_vac[1] + step
                ok = 0
                move_y = 1
            else:
                move_y = 0
        
            #if z_min < 1.5: #box_z / 3.:
            if (len(self.atomsInLayer(atomsTemp, 2, 0., 0.25)) != 8) or (len(self.atomsInLayer(atomsTemp, 2, 6.25, 0.25)) != 8):
                trans_vac[2] = trans_vac[2] + step
                ok = 0
                move_z = 1
            else:
                move_z = 0
        
#            print trans_vac
            if ok:
                
                #                atomsTemp.translate([step * move_x, step * move_y, step * move_z])
                #atomsTemp.set_scaled_positions(atomsTemp.get_scaled_positions())
                posTrans = atomsTemp.get_positions()[self.missingNeighbors]

                x_vac = np.average(posTrans[:,0])
                y_vac = np.average(posTrans[:,1])
                z_vac = np.average(posTrans[:,2])
                pos_vac = np.array([x_vac, y_vac, z_vac])
#
#                a1 = Atom('O', pos_vac)
#                view(atomsTemp + a1)
#                
#                a2 = Atom('S', pos_vac - trans_vac)
#                view(atomsTempOrg + a2)
#                
#                raw_input()
                break
    
        if not ok:
            pos_vac = trans_vac
            print 'problems'
            return None
        return pos_vac - trans_vac
            
    def particleToParticleVec(self, p1, p2):
        D = p2 -p1
        Dr = la.solve(self.atoms._cell.T, D)
        return np.dot(Dr - np.round(Dr) * self.atoms._pbc, self.atoms._cell)
        
    def display(self):
        #a = Atom('O', self.calcMeanLocation())
        a = Atom('O', self.calcMeanLocation_new())
        view(self.atoms + a)

    def atomsInLayer(self, atoms, index, value, width):
        pos = atoms.get_positions()
        below = (pos[:,index] < (value + width))
        above = (pos[:,index] > (value - width))

#        print index
#        print pos[:,index]
#        print below & above
#        print sum (below & above)
#        raw_input()
        
        return [layer[0] for layer in enumerate(below & above) if layer[1]]
        
#        return sum( below & above )
