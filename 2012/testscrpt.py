import _path
import numpy as np
from ChangingNeighbors import *
from SimulationData import *
from enthought.mayavi import mlab

path = '/home/kjartan/Temp/50Zr_50Cu'
sData = SimulationData(path)
cutoffs = {1:{1:3.3,2:3.7},2:{1:3.7,2:3.9}}


def plotThisThing(state, proc):
    react = sData.getInitReactant(state)
    prod = sData.getInitReactant(proc)
    #react, prod = sData.getReacProd(state,proc)
    cn = ChangingNeighbors(react, prod, cutoffs)


    ln = cn.nrValues(cn.lostNeighbors)
    indxs = np.where(ln > 0)
    ln = ln[indxs]
    pos = react.get_positions()[indxs]
    t = react.get_atomic_numbers()[indxs]
    pts = mlab.quiver3d(pos[:,0], pos[:,1], pos[:,2],t,t,t,scalars=ln, mode='sphere')

    pts.glyph.color_mode = 'color_by_scalar'
    pts.glyph.glyph_source.glyph_source.center = [0,0,0]

    indx2 = np.setdiff1d(np.array(range(prod.get_number_of_atoms())),indxs[0])
    p = react.get_positions()[indx2]
    t2 = react.get_atomic_numbers()[indx2]
    pt = mlab.quiver3d(p[:,0],p[:,1],p[:,2],t2,t2,t2, opacity=.06, mode='sphere')
    
    mlab.figure()
    
    ln = cn.nrValues(cn.newNeighbors)
    indxs = np.where(ln > 0)
    ln = ln[indxs]
    pos = react.get_positions()[indxs]
    t = react.get_atomic_numbers()[indxs]
    pts = mlab.quiver3d(pos[:,0], pos[:,1], pos[:,2],t,t,t,scalars=ln, mode='sphere')

    pts.glyph.color_mode = 'color_by_scalar'
    pts.glyph.glyph_source.glyph_source.center = [0,0,0]

    indx2 = np.setdiff1d(np.array(range(prod.get_number_of_atoms())),indxs[0])
    p = react.get_positions()[indx2]
    t2 = react.get_atomic_numbers()[indx2]
    pt = mlab.quiver3d(p[:,0],p[:,1],p[:,2],t2,t2,t2, opacity=.06, mode='sphere')
    
    mlab.figure()
    
    ln = cn.nrValues(cn.lostNeighbors) + cn.nrValues(cn.newNeighbors)
    indxs = np.where(ln > 0)
    ln = ln[indxs]
    pos = react.get_positions()[indxs]
    t = react.get_atomic_numbers()[indxs]
    pts = mlab.quiver3d(pos[:,0], pos[:,1], pos[:,2],t,t,t,scalars=ln, mode='sphere')

    pts.glyph.color_mode = 'color_by_scalar'
    pts.glyph.glyph_source.glyph_source.center = [0,0,0]

    indx2 = np.setdiff1d(np.array(range(prod.get_number_of_atoms())),indxs[0])
    p = react.get_positions()[indx2]
    t2 = react.get_atomic_numbers()[indx2]
    pt = mlab.quiver3d(p[:,0],p[:,1],p[:,2],t2,t2,t2, opacity=.06, mode='sphere')
    
