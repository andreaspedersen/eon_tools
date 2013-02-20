import multiprocessing
from multiprocessing import Pool
import _path
import numpy as np
from SimulationData import *
from FindNeighbors import *
import sys

#path = '/data/users/kgm/raw_data/30Zr_70Cu'
path = '/home/andreas/simulations/eon-cuzr/temp'

sData = SimulationData(path)

def mainfunc(i):
    f = sData.getInitialReac(i)
    if f:
        FindNeighbors(f)
    print 'State %i is finished' % (i,)

if __name__ == '__main__':
    nrCPU = multiprocessing.cpu_count()
    print 'The number of processes will be %i' % (nrCPU,)

    pool = Pool(processes=nrCPU)
    pool.map(mainfunc, xrange(sData.nrStates))
