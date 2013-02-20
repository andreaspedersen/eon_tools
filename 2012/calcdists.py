import _path
import numpy as np
from SimulationData import *
from FindNeighbors import *
import sys

if len(sys.argv) < 3:
    exit()

#path = '/data/users/kgm/raw_data/30Zr_70Cu'
path = '/home/kjartan/Temp/50Zr_50Cu'

sData = SimulationData(path)

for i in xrange(int(sys.argv[1]), int(sys.argv[2])):
    f = sData.getInitialReac(i)
    if f:
        FindNeighbors(f)
    else:
        print 'File for state %i not found.' % (i,)
print 'Finished states %i to %i' %(int(sys.argv[1]), int(sys.argv[2]))
