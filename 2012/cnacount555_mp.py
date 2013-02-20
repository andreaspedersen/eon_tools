import multiprocessing
from multiprocessing import Pool
import _path
import numpy as np
from SimulationData import *
from FindNeighbors import *
import sys
import cna

#path = '/data/users/kgm/raw_data/30Zr_70Cu'
path = '/home/andreas/simulations/eon-cuzr'
#path = '/home/kjartan/Temp/50Zr_50Cu'

sData = SimulationData(path)

def mainfunc(i):
    f = sData.getInitialReac(i)
    if f:
        fn = FindNeighbors(f)
        fn.setCutOffs({1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})
        cn = cna.full_cna(fn)
        f2 = f[:-4] + '_cna.txt'
        np.savetxt(f2, cn)
        return cna.get_unique_cna(cn).get((5,5,5),0)

if __name__ == '__main__':
    nrCPU = multiprocessing.cpu_count()
    #print 'The number of processes will be %i and the number of states are %i' % (nrCPU, nrStates)

    pool = Pool(processes=nrCPU)
    res = pool.map(mainfunc, xrange(sData.nrStates))
    
    np.savetxt(path + '/555s.txt', np.array(res,dtype='int'))
