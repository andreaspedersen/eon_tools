# use script in cmapdata in thor to generate the needed files

import _path
import numpy as np
import matplotlib.pyplot as plt
import os
import AtomsTools as at
import Converter as cc
from SimulationData import *

path = "/home/kjartan/Temp/out3"

#energy of SP versus X, Y (i and j) of the atom that moved the most
pos = np.loadtxt(os.path.join(path, "positions.txt"))
sadEn = np.loadtxt(os.path.join(path, "saddEnergy.txt"))
cell = 27.688400

pos[pos>cell] = 0

# i and j represent the index of the axis
def doColormap(i,k,di,dk):
    d = np.array([di,dk])
    x = int(np.ceil(cell/di))
    y = int(np.ceil(cell/dk))
    cm = np.ones((x,y)) * max(sadEn)
    
    for j in xrange(sadEn.size):
        c = map(int, np.floor(pos[j,(i,k)]/d))
        #if cm[c[0],c[1]] != 0:
        #    print 'Over writing at pos ', pos[j]
        cm[c[0],c[1]] = sadEn[j]
    
    plt.pcolor(cm)
    return cm
        
    
