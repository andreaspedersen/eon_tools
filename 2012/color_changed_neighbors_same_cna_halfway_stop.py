import _path
import numpy as np
import os
from find_neighbors import *
from changing_neighbors import *
import coloring_atoms as ca
import converter as cc
from ase import view

path = '/home/kjartan/Temp/chneigh'

chn1 = ChangingNeighbors(os.path.join(path, '0', 'reactant.con'),
                         os.path.join(path, '2387', 'reactant.con'),
                         {1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})
                         
chn2 = ChangingNeighbors(os.path.join(path, '0', 'reactant.con'),
                         os.path.join(path, '1635', 'reactant.con'),
                         {1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})

chn3 = ChangingNeighbors(os.path.join(path, '1635', 'reactant.con'),
                         os.path.join(path, '2387', 'reactant.con'),
                         {1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})
                         
reac1 = chn1.neigh1.atoms

t1 = chn1.nrValues(chn1.lostNeighbors)
t2 = chn1.nrValues(chn2.lostNeighbors)
t3 = chn1.nrValues(chn3.lostNeighbors)

indx1 = (t1 > 0)
indx2 = (t2 > 0)
indx3 = (t3 > 0)

aseAtoms = cc.eon_2_ase(reac1)

aseAtoms = ca.color(aseAtoms,indexes=indx1,color='red', makeCopy=False)
aseAtoms = ca.color(aseAtoms,indexes=indx2,color='green', makeCopy=False)
aseAtoms = ca.color(aseAtoms,indexes=indx3,color='blue', makeCopy=False)
aseAtoms = ca.color(aseAtoms,indexes=np.logical_and(indx2,indx3),color='yellow', makeCopy=False)

view(aseAtoms)

