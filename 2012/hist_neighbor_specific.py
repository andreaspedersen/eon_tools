import _path
from SimulationData import *
from FindNeighbors import *
import numpy as np
import matplotlib.pyplot as plt
import os

cutoffs = np.array([[3.9,3.7],[3.7,3.3]])
outPath = '/home/kjartan/Temp/out'

sData = SimulationData('/home/kjartan/Temp/50Zr_50Cu')
fn = FindNeighbors(sData.getInitReactant(0))
#x = fn.rdf(9.,9)

"""
for i, atoms in enumerate(sData.getAllInitReactant()):
    fn = FindNeighbors(atoms)
    data = []
    
    for aT, nT in ((1,1),(1,2),(2,1),(2,2)):
        data.append(fn.countNeighborsByType(aT,nT,cutoffs[aT-1,nT-1]))
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ba1 = ax.bar(range(data[0].size),data[0],width=1.,color='red',alpha=.3)
    ba2 = ax.bar(range(data[1].size),data[1],width=1.,alpha=.3)
    fig.savefig(os.path.join(outPath, "figure%i_1.png" % (i,)))
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ba1 = ax.bar(range(data[2].size),data[2],width=1.,color='red',alpha=.3)
    ba2 = ax.bar(range(data[3].size),data[3],width=1.,alpha=.3)
    fig.savefig(os.path.join(outPath, "figure%i_2.png" % (i,)))
"""
