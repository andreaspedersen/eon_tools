#!/usr/bin/python
import numpy
import sys
import math

class VerletCell():
    def __init__(self,X1,X2,Y1,Y2,Z1,Z2):
        self._list = []
        assert(X1<X2 and Y1<Y2 and Z1<Z2)
        self._x1 = X1
        self._x2 = X2
        self._y1 = Y1
        self._y2 = Y2
        self._z1 = Z1
        self._z2 = Z2

    def clean(self):
        self._list=[]

    def append(self,entry):
        self._list.append(entry)

    def listEntries(self):
        for entry in self._list:
            yield entry

    def numEntries(self):
        return len(self._list)

    def remove(self,entry):
        lenbefore = self.numEntries()
        self._list = [ i for i in self._list if i !=entry ]
        if self.numEntries() != lenbefore-1:
            sys.stderr.write("WARNING: remove() operation on verletCell removed %d entries\n" % (lenbefore - self.numEntries() ) )
       
    def __iter__(self):
        for entry in self._list:
            yield entry


class VerletList():
    def __init__(self):
        self._cellList = []
        self._neighbourList = []
        self._box = numpy.zeros(3)
        self._numCells = None
        self._ncell_v = numpy.zeros(3,dtype=int)
        self._cellLength = numpy.zeros(3)
        self._entries= { }
        self._shift = numpy.zeros(3)

    def initialize(self,box,rc,shift = None):
        '''
            Create cell list for a box with cuttof rc
        '''
        self._box = numpy.diag(box)        
        for i in range(3):
            print "box", self._box[i]
            self._ncell_v[i] = int(math.floor(self._box[i]/rc))
            if self._ncell_v[i] < 3:
                self._ncell_v[i]=3
        
        
        self._numCells = self._ncell_v[0]*self._ncell_v[1]*self._ncell_v[2]
        
        self._cellLength = self._box / self._ncell_v
        nx = self._ncell_v[0]
        ny = self._ncell_v[1]
        nz = self._ncell_v[2]

        #create cell
        icell = 0        
        for iz in range(self._ncell_v[2]):
            for iy in range(self._ncell_v[1]):
                for ix in range(self._ncell_v[0]):
                    cell = VerletCell(  ix*self._cellLength[0]+self._shift[0],(ix+1)*self._cellLength[0]+self._shift[0],
                                        iy*self._cellLength[1]+self._shift[1],(iy+1)*self._cellLength[1]+self._shift[1],
                                        iz*self._cellLength[2]+self._shift[2],(iz+1)*self._cellLength[2]+self._shift[2] 
                                     )
                    nblist = []
                    #calculate nearest neighbours (6)
                    #left
                    nbid = icell - 1 if (ix != 0) else icell+ (nx-1)           
                    nblist.append(nbid)
                    #right                 
                    nbid = icell + 1 if (ix != nx-1) else icell - (nx-1) 
                    nblist.append(nbid)
                    #front                    
                    nbid = icell - nx if (iy != 0) else icell + (ny-1)*nx
                    nblist.append(nbid)
                    #back                    
                    nbid = icell + nx if (iy!=ny-1) else icell - (ny-1)*nx 
                    nblist.append(nbid)
                    #bottom
                    nbid = icell - nx*ny if (iz!=0) else icell + (nz-1)*nx*ny
                    nblist.append(nbid)
                    #top
                    nbid = icell + nx*ny if (iz!=nz-1) else icell - (nz-1)*nx*ny
                    nblist.append(nbid)

                    self._cellList.append(cell)
                    self._neighbourList.append(nblist)

                    icell+=1

        #calculate other neighbours
        for icell in range(self._numCells):
            #print self._neighbourList[icell]
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][0] ][2] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][0] ][3] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][1] ][2] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][1] ][3] )

            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][0] ][4] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][0] ][5] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][1] ][4] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][1] ][5] )

            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][2] ][4] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][2] ][5] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][3] ][4] )
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList[icell][3] ][5] )

            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][0] ][2] ][4])
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][0] ][3] ][4])
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][1] ][2] ][4])
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][1] ][3] ][4])

            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][0] ][2] ][5])
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][0] ][3] ][5])
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][1] ][2] ][5])
            self._neighbourList[icell].append( self._neighbourList [ self._neighbourList [ self._neighbourList[icell][1] ][3] ][5])


    def clean(self):
        for cell in self._cellList:
            cell.clean()
        self._entries = { }    
            

    def neighbours(self,entry):
        try:
            cellid = self._entries[entry]
        except KeyError:
            print "Cannot find %d in verlet list to give neighbours. Entry not present." % entry
            return         
            
        for entry in self._cellList[cellid]:
            yield entry
        for neighbourcellid in self._neighbourList[cellid]:
            for entry in self._cellList[neighbourcellid]:
                yield entry

    def cell(self,entry):
        try:
            cellid = self._entries[entry]
        except KeyError:
            print "ERROR"
            return
        return self._cellList[cellid]                    

    def remove(self,entry):
        try:
            cellid = self._entries[entry]
        except KeyError:
            print "Cannot remove %d from verlet list. Entry not present." % entry
            return
            
        self._cellList[cellid].remove(entry)
        self._entries.pop(entry)

    def append(self,entry,coords):
        if entry in self._entries.keys():
            print "ERROR: %d already in verlet list" % entry
            return
            
        #calculate verlet cell:
        ix = 0
        while coords[0] >= (ix+1)*self._cellLength[0]:     
            ix+=1                
        iy = 0
        while coords[1] >= (iy+1)*self._cellLength[1]:     
            iy+=1   
        iz = 0
        while coords[2] >= (iz+1)*self._cellLength[2]:     
            iz+=1
        
        cellid = ix + self._ncell_v[0]*iy + self._ncell_v[0]*self._ncell_v[1]*iz
        self._cellList[cellid].append(entry)
        self._entries[entry]=cellid
                       
