# not useful

import _path
from SimulationData import *
from FindNeighbors import *
import numpy as np
import matplotlib.pyplot as plt
import os

sData = SimulationData('/home/kjartan/Temp/50Zr_50Cu')
fn = FindNeighbors(sData.getInitReactant(0), {1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})
result = zeros(12,dtype='int')
for i in xrange(12):
    neigh = np.hstack((fn.getNeighbors(i,1),fn.getNeighbors(i,2)))
    nneigh = zeros(0, dtype='int')
    print 'I is here!'

    for neighbor in neigh:
        for j in xrange(1,3):
            nneigh = np.unique(np.hstack((nneigh, fn.getNeighbors(neighbor, j))))
    result[i] = nneigh.size
