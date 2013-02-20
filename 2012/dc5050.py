# probably
# find the neighbors for the lowest product state
import _path
import numpy as np
import matplotlib.pyplot as plt
import os
import AtomsTools as at
import Converter as cc
from SimulationData import *
from FindNeighbors import *

path = '/home/kjartan/Temp/2/50Cu50Zr/981/reactant'
#sData = SimulationData(path)
#energy, states, process = sData.getLowestProductEnergyOfAll(10)

fn = FindNeighbors(path)
fn.setCutOffs({1:{1:3.3,2:3.7},2:{1:3.7,2:3.9}})

#diff = np.zeros(energy.size)
#for i in xrange(energy.size):
#    diff[i] = energy[i] - sData.getSaddleEnergys(states[i])[i]
