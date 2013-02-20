import _path
import numpy as np
import os
from find_neighbors import *
from changing_neighbors import *

path = '/home/kjartan/Temp/chneigh'

chn = ChangingNeighbors(os.path.join(path, '0', 'reactant.con'),
                        os.path.join(path, '2387', 'reactant.con'),
                        {1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})
reac1 = chn.neigh1.atoms

t = chn.nrValues(chn.lostNeighbors)
an = reac1.get_atomic_numbers()

indx1 = (t > 0)
indx2 = (t == 0)

q1 = reac1.r[indx1]
q2 = reac1.r[indx2]

