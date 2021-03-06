"""
ATTENTION: This file is here for legacy reasons only, the code has been moved
to eon_tools under the eon_read file, use that.
"""


"""
Author: Kjartan Myrdal
Date: June 2012
License: Undefined, if you have this you're probably allowed to use it.
"""

import os, sys
import numpy as np
import numpy.linalg as la
import ConfigParser

"""
Takes a path to data created by eOn AKMC calculation and makes the data simply
accesible through functioncalls.
"""
class SimulationData(object):

    def __init__(self, path='.'):
        self.path = path
        self.nrStates = self._findNumStates()
        self.nrProcesses = self._findNumProc()
        self.nrTabulatedSearches = self._findTabSearches()
        self.totalNrProcesses = self._findTotProc()
        
    def _findNumStates(self):
        fName = 'state_table'
        path = os.path.join(self.path, 'states')
        try:
            res = len(np.loadtxt(os.path.join(path, fName), usecols=(0,)))
        except IOError:
            print 'file state_table is empty. Assuming single state.'
            res = 1
        return res
        
    def _findNumProc(self):
        procs = np.empty(self.nrStates, dtype='int')
        fName = 'processtable'
        for state in xrange(self.nrStates):
            path = os.path.join(self.path, 'states', str(state))
            out = np.loadtxt(os.path.join(path, fName),
                             skiprows=1, usecols=(0,))
            procs[state] = out.size
        return procs

    def _findTabSearches(self):
        procs = np.empty(self.nrStates, dtype='int')
        fName = 'processtable'
        for state in xrange(self.nrStates):
            path = os.path.join(self.path, 'states', str(state))
            out = np.loadtxt(os.path.join(path, fName),
                             skiprows=1, usecols=(8,))
                             #print state
                             #print np.size(out)
            #procs[state] = out.size
            if np.size(out) > 1:
                procs[state] = int(sum(out))
            else:
                procs[state] = int(out)
            #            print out
            #print int(sum(out))
            #raw_input()
        return procs

    def _findTotProc(self):
        result = 0
        for key in xrange(self.nrStates):
            result += self.nrProcesses[key]
        return result
     
    def _getProcessTable(self, state):
        if state < self.nrStates:
            path = os.path.join(self.path, 'states', str(state), 'processtable')
            return path
    
    def _getOrderBy(self, f, index):
        """
        Takes column in a file (proctable) and returns sorted indexed
        """
        arr = np.loadtxt(f,skiprows=1,usecols=(index,))
        return np.argsort(arr)
    
    def _stateProcToIndex(self, state, proc):
        return self.nrProcesses[:state].sum() + proc
    
    def _indexToStateProc(self, index):
        states = 0
        if index < 0:
            raise IndexError("Negative indexes not allowed.")
        for state in xrange(self.nrStates):
            
            if index < states + self.nrProcesses[state]:
                proc = index - states
                return state, proc
            states += self.nrProcesses[state]
        raise IndexError("Index %i is out of range. Highest value is %i" %(index, states))
    
    
    ## functions addressing a single specified state
    
    def forStateIdGetProcessIndexWithLowestSaddleEnergy(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = self._getOrderBy(f, 1)
        f.close()
        return res
    
    def forStateIdGetProcessIndexWithLowestProductEnergy(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = self._getOrderBy(f, 4)
        f.close()
        return res
    
    def forStateIdGetReactantEnergy(self, state):
        filePath = os.path.join(self.path, 'states', str(state), 'info')
        f = open(filePath, 'r')
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(f)
        f.close()
        return cfg.getfloat('MetaData', 'reactant energy')

    #----------

    def forStateIdGetSaddleEnergies(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(1,))
        f.close()
        if np.size(res) == 1:
            res = np.array([res])
        return res
    
    def forStateIdGetProductEnergies(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(4,))
        f.close()
        if np.size(res) == 1:
            res = np.array([res])
        return res

    def forStateIdGetBarriers(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(6,))
        f.close()
        if np.size(res) == 1:
            res = np.array([res])
        return res

    def forStateIdGetPrefactors(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(2,))
        f.close()
        if np.size(res) == 1:
            res = np.array([res])
        return res

    def forStateIdGetRates(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(7,))
        f.close()
        if np.size(res) == 1:
            res = np.array([res])
        return res

    # returns a list with the product ids for the different processes
    # (-1 = not traversed)
    def forStateIdGetProductsAccepted(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(3,))
        f.close()
        if np.size(res) == 1:
            res = np.array([res])
        return res

    #----------


    def forStateIdGetAcceptedSaddleEnergies(self, state):
        res = self.forStateIdGetSaddleEnergies(state)
        picked = self.forStateIdGetProductsAccepted(state) > -1
        return res[np.where(picked)[0]]
    
    def forStateIdGetAcceptedProductEnergies(self, state):
        res = self.forStateIdGetProductEnergies(state)
        picked = self.forStateIdGetProductsAccepted(state) > -1
        return res[np.where(picked)[0]]
    
    def forStateIdGetAcceptedBarriers(self, state):
        res = self.forStateIdGetBarriers(state)
        picked = self.forStateIdGetProductsAccepted(state) > -1
        return res[np.where(picked)[0]]
    
    def forStateIdGetAcceptedPrefactors(self, state):
        res = self.forStateIdGetPrefactors(state)
        picked = self.forStateIdGetProductsAccepted(state) > -1
        return res[np.where(picked)[0]]

    # functions addressing all states visited
    
    def fromAllStatesGetProductEnergies(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetProductEnergies(i)))
        return np.array(l,'f')

    def fromAllStatesGetSaddleEnergies(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetSaddleEnergies(i)))
        return np.array(l,'f')

    def fromAllStatesGetBarriers(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetBarriers(i)))
        return np.array(l,'f')

    def fromAllStatesGetPrefactors(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetPrefactors(i)))
        return np.array(l,'f')

    #----------
    
    def fromAllStatesGetAcceptedProductEnergies(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetAcceptedProductEnergies(i)))
        return np.array(l,'f')
    
    def fromAllStatesGetAcceptedSaddleEnergies(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetAcceptedSaddleEnergies(i)))
        return np.array(l,'f')
    
    def fromAllStatesGetAcceptedBarriers(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetAcceptedBarriers(i)))
        return np.array(l,'f')
    
    def fromAllStatesGetAcceptedPrefactors(self):
        l = []
        for i in range(self.nrStates):
            l.extend(list(self.forStateIdGetAcceptedPrefactors(i)))
        return np.array(l,'f')

#----------

    def fromAllStatesGetLowestEnergyProducts(self, nrProc=-1):
        """
        returns a slice of the lowest products up until the index 'nrProc'
        """
        energy = np.zeros(0)
        states = np.zeros(0, dtype='int')
        process = np.zeros(0, dtype='int')
        for state in xrange(self.nrStates):
            f = open(self._getProcessTable(state), 'r')
            arr = np.loadtxt(f,skiprows=1,usecols=(4,))
            f.close()
            energy = np.hstack((energy,arr))
            arg = np.argsort(energy)
            states = np.hstack((states,np.ones(arr.size,dtype='int')*state))
            process = np.hstack((process,
                                 np.array(range(arr.size),dtype='int')))
            energy = energy[arg][:nrProc]
            states = states[arg][:nrProc]
            process = process[arg][:nrProc]
        return energy, states, process

    
    # functions to generate paths
    
#    def getReac(self, state, proc):
    def getPathProcessReactant(self, state, proc):
        """
        get path for specific reactant for specific state
        """
        if state < self.nrStates and proc < self.nrProcesses[state]:
            fName = 'reactant_%i.con' %(proc,)
            path = os.path.join(self.path, 'states', str(state), 'procdata')
            return os.path.join(path, fName)
        print "State or process doesn't exist"
    
#    def getProd(self, state, proc):
    def getPathProcessProduct(self, state, proc):
        """
        get path for specific product for specific state
        """
        if state < self.nrStates and proc < self.nrProcesses[state]:
            fName = 'product_%i.con' %(proc,)
            path = os.path.join(self.path, 'states', str(state), 'procdata')
            return os.path.join(path, fName)
        print "State or process doesn't exist"
    
#    def getReacProd(self, state, proc):
    def getPathsProcessReactantAndProduct(self, state, proc):
        """
        get paths for specific process for specific state
        """
        return self.getPathProcessReactant(state,proc), self.getPathProcessProduct(state,proc)
    
#    def getAllReacProd(self):
    def getPathsAllProcessReactantAndProduct(self):
        """
        get paths for all processes in all proctables
        """
        for state in xrange(self.nrStates):
            for proc in xrange(self.nrProcesses[state]):
#                yield self.getReacProd(state, proc), (state, proc)
                yield self.getPathsProcessReactantAndProduct(state, proc), (state, proc)
            
#   def getInitialReac(self, state):
    def getPathStateReactant(self, state):
        if state < self.nrStates:
            fName = 'reactant.con'
            path = os.path.join(self.path, 'states', str(state))
            return os.path.join(path, fName)
        print "State doesn't exist"
                

    
    
        
