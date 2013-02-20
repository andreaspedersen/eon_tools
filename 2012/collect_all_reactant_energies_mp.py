import multiprocessing
from multiprocessing import Pool
import _path
import numpy as np
from SimulationData import *
import sys, os

#path = '/data/users/kgm/raw_data/30Zr_70Cu'
#path = '/home/andreas/simulations/eon-cuzr'
path = '/home/kjartan/Temp/50Zr_50Cu'

sData = SimulationData(path)

def mainfunc(i):
    f = sData.getReactantEnergy(i)
    if f:
        return f
        
if __name__ == '__main__':
    nrCPU = multiprocessing.cpu_count()
    nrCPU = 4
    #print 'The number of processes will be %i and the number of states are %i' % (nrCPU, nrStates)

    pool = Pool(processes=nrCPU)
    res = pool.map(mainfunc, xrange(sData.nrStates))
    
    np.savetxt(os.path.join(path, 'reacEns.txt'), np.array(res))
