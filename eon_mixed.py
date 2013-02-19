import eon_read

#from pylab import *
#from eon_analyze import *
#from eon_extract import *

def project_to_plane(traj,index_0,index_1):
    x = []
    y = []
    for state in traj:
        x.append(state[index_0])
        y.append(state[index_1])
    return x,y

def make_trac_in_plane(diff,index_0,index_1):
    x = [0]
    y = [0]
    i = 0
    for state in diff:
        x.append(x[i] + state[index_0])
        y.append(y[i] + state[index_1])
        i = i+1
    return x,y

