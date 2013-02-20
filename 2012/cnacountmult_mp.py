import multiprocessing
from multiprocessing import Pool
import _path
import numpy as np
from SimulationData import *
from FindNeighbors import *
import sys, os
import cna

#path = '/data/users/kgm/raw_data/30Zr_70Cu'
#path = '/home/andreas/simulations/eon-cuzr'
path = '/home/kjartan/Temp/50Zr_50Cu'

sData = SimulationData(path)

def mainfunc(i):
    f = sData.getInitialReac(i)
    if f:
        f2 = f[:-4] + '_cna.txt'
        cn = np.loadtxt(f2)
        ucn = cna.get_unique_cna(cn)
        return [ucn.get(k,0) for k in ((4,2,1),(4,1,1),(6,6,6))]

if __name__ == '__main__':
    nrCPU = multiprocessing.cpu_count()
    nrCPU = 2
    #print 'The number of processes will be %i and the number of states are %i' % (nrCPU, nrStates)

    pool = Pool(processes=nrCPU)
    res = np.array(pool.map(mainfunc, xrange(sData.nrStates)))
    
    np.savetxt(os.path.join(path, '421s.txt'), res[:,0],fmt='%i')
    np.savetxt(os.path.join(path, '411s.txt'), res[:,1],fmt='%i')
    np.savetxt(os.path.join(path, '666s.txt'), res[:,2],fmt='%i')
