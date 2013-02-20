# colors the atoms with same local environment (cna) blue
# colors the atoms with changed neighbors red
# grey atoms same neighbors different cna


import _path
import numpy as np
import os
from find_neighbors import *
from changing_neighbors import *
import coloring_atoms as ca
import converter as cc
from ase import view
import cna

path = '/home/kjartan/Temp/chneigh'

chn1 = ChangingNeighbors(os.path.join(path, '0', 'reactant.con'),
                         os.path.join(path, '984', 'reactant.con'),
                         {1:{1:3.9,2:3.7},2:{1:3.7,2:3.3}})

reac1 = chn1.neigh1.atoms

t1 = chn1.nrValues(chn1.lostNeighbors)
t2 = chn1.nrValues(chn1.newNeighbors)

# indecies of atoms where neighbors are changing
indx1 = (np.logical_or(t1 > 0,t2 > 0))
indx1i = np.arange(indx1.size)[indx1]

aseAtoms = cc.eon_2_ase(reac1)

aseAtoms = ca.color(aseAtoms,indexes=indx1,color='red', makeCopy=False)

cn1 = cna.full_cna(chn1.neigh1)
cn2 = cna.full_cna(chn1.neigh2)
#cn1 = np.loadtxt('/home/kjartan/Temp/chneigh/0/reactant_cna.txt',dtype='int')
#cn2 = np.loadtxt('/home/kjartan/Temp/chneigh/984/reactant_cna.txt',dtype='int')

cn1i = cna.get_cna_by_indexes(indx1i,cn1)
cn2i = cna.get_cna_by_indexes(indx1i,cn2)

u1 = cna.get_unique_cna(cn1i)
u2 = cna.get_unique_cna(cn2i)

# indecies of atoms with same cna
indx2 = cna.cna_atoms_no_change(cn1, cn2)

aseAtoms = ca.color(aseAtoms,indexes=indx2,color='blue', makeCopy=False)

view(aseAtoms)

