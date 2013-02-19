import numpy
import sys
import verletlist

def length_angle_to_box(boxlengths, angles):
    box = numpy.zeros( (3,3) )
    angles *= numpy.pi/180.0
    box[0][0] = 1.0
    box[1][0] = numpy.cos(angles[0])
    box[1][1] = numpy.sin(angles[0])
    box[2][0] = numpy.cos(angles[1])
    box[2][1] = (numpy.cos(angles[2])-box[1][0]*box[2][0])/box[1][1]
    box[2][2] = numpy.sqrt(1.0-box[2][0]**2-box[2][1]**2)
    box[0,:]*=boxlengths[0]
    box[1,:]*=boxlengths[1]
    box[2,:]*=boxlengths[2]
    return box


class Atoms(object):
    """ The Atoms class. """

    def __init__(self, atoms):
        if type(atoms)  == str:
            self.read_con(atoms)
        else:
            self.r = numpy.zeros((atoms,3))
            self.free = numpy.ones(atoms)
            self.box = numpy.zeros((3,3))
            self.names = ['']*atoms
            self.mass = numpy.zeros(atoms)
            self.num_atoms = atoms
            self.verletList = None
            self.verletCutoff = None

    def __len__(self):
        '''
        Returns the number of atoms in the object'''
        #return len(self.r)
        return self.num_atoms
    
    def get_atomic_numbers(self):
        return numpy.vectorize(lambda x: elements[x]['number'])(self.names)
    
    def set_atomic_numbers(self, numbers):
        self.names = numpy.vectorize(lambda x: elements[x]['symbol'])(numbers)

    # functions defined above
    atomic_numbers = property(get_atomic_numbers,set_atomic_numbers)
        
    def copy(self):
        p = Atoms(len(self))
        p.r = self.r.copy()
        p.free = self.free.copy()
        p.box = self.box.copy()
        p.names = self.names[:]
        p.mass = self.mass
        p.num_atoms = self.num_atoms
        return p

    def free_r(self):
        nfree = sum(self.free)
        temp = numpy.zeros((nfree, 3))
        index = 0
        for i in range(len(self.r)):
            if self.free[i]:
                temp[index] = self.r[i]
                index += 1
        return temp

    def appendvec(self,atomicnr,r,free=1):
        self.append(atomicnr,r[0],r[1],r[2],free)
    
    def append(self,atomicnr,x,y,z,free=1):
        if atomicnr.__class__==str:
            #find the atomic nr.
            ielement = 0
            for key in elements:
                if elements[key]['symbol']==atomicnr:
                    break
                ielement+=1
            atomicnr=ielement
        
        if atomicnr==-1 or atomicnr=="X":
            return
        self.num_atoms +=1
        newatom=numpy.array( [ float(x),float(y),float(z) ] )
        self.r=numpy.vstack((self.r,newatom))
        self.free=numpy.concatenate( (self.free, [ int(free) ] ) )
        self.mass = numpy.concatenate ( ( self.mass,[ elements[int(atomicnr)]['mass'] ] ) )
        self.names.append( elements[int(atomicnr)]['symbol'] )
        return
    
    def elements(self):
        """
        Determine dictionary with the number of different elements
        """        
        d=dict()
        for a_type in numpy.unique(self.names):
            d[a_type] = self.names[self.names == a_type]
        return d
    """
    def elements(self):
        '''
        Determine dictionary with the number of different elements
        '''
        d=dict()
        for iatom in range(len(self)):
            if self.names[iatom] not in d:
                d[self.names[iatom] ] = 0
        
            d[self.names[iatom] ]+=1
        return d
    """

    def create_verletList(self,cutoff):
        if self.verletList != None:
            if self.verletCutoff == cutoff:
                self.verletList.clean()
                return
            del self.verletList
        
        shift = self.r[0]
        for iatom in range(self.num_atoms):
            r = self.r[iatom]
            for i in range(3):
                if r[i] < shift[i]:
                    shift[i] = r[i]
        
        self.verletList = verletlist.VerletList()
        self.verletList.initialize(self.box,cutoff,shift)
        self.verletCutoff = cutoff
    
    def get_vector_points(self, p1, p2):
        ibox = numpy.linalg.inv(self.box)
        return pbc(p2 - p1, self.box, ibox)
    
    def get_distance_points(self, p1, p2):
        return numpy.linalg.norm(self.get_vector_points(p1, p2))
    
    def get_vector_atoms(self, i, j):
        return self.get_vector_points(self.r[i], self.r[j])
    
    def get_distance_atoms(self, i, j):
        return self.get_distance_points(self.r[i], self.r[j])

#################################################################################################
#  Storage related stuff
#################################################################################################
    def read_con(self,filein):
        '''
        Load a con file
            filein: may be either a filename or a file-like object
        '''
        if hasattr(filein, 'readline'):
            con = filein
        else:
            con = open(filein, 'r')
        
        con.readline() #Line 1: comment
        con.readline() #Line 2: comment
        
        boxlengths = numpy.fromstring(con.readline(), sep=' ') #Line 3: Box lengths
        dim = boxlengths.size
        boxangles = numpy.fromstring(con.readline(),count=dim, sep=' ')
        if not ((abs(boxangles - 90) < 1e-16).all()):
            boxtemp = length_angle_to_box(boxlengths,boxangles)
        else:
            boxtemp = numpy.diag(boxlengths)
        
        con.readline() #Line 5: comment
        con.readline() #Line 6: comment
        
        num_types = int(con.readline().strip())
        num_each_type = numpy.fromstring(con.readline(),dtype='int', sep=' ')
        mass_of_type = numpy.fromstring(con.readline(), sep=' ')
        num_atoms =  num_each_type.sum()
        self.num_atoms = num_atoms
        self.box = boxtemp
        self.r = numpy.empty((num_atoms,dim))
        self.free = numpy.empty(num_atoms, dtype=int)
        self.names = numpy.empty(num_atoms, dtype='|S2')
        self.mass = numpy.empty(num_atoms)
        
        j_sum = 0
        for i in xrange(num_types):
            i_sum = j_sum
            j_sum += num_each_type[i]
            name = con.readline().strip()
            self.names[i_sum:j_sum] = name
            self.mass[i_sum:j_sum] = mass_of_type[i]
            con.readline()
            for j in xrange(i_sum, j_sum):
                line = numpy.fromstring(con.readline(), count=dim+1, sep=' ')
                self.r[j] = line[:dim]
                #Con file stores bool of fixed atoms this file uses bool of free
                self.free[j] = abs(int(line[-1] - 1))
    """
        def read_con(self,filein):
        
        '''
        Load a con file
        filein: may be either a filename or a file-like object
        '''
        if hasattr(filein, 'readline'):
        con = filein
        else:
        con = open(filein, 'r')
        con.readline() #Line 1: comment
        con.readline() #Line 2: comment
        # determine how many dimensions
        tmp = numpy.array(con.readline().split()) #Line 3: Box lengths
        for i in range(len(tmp)):
            dim=i+1
            try: float(tmp[i])
            except:
                dim=i
                break
        #handle the box
        boxlengths=numpy.zeros(dim)
        for i in range(dim):
            boxlengths[i]=float(tmp[i])
        boxangles=numpy.array([ float(f) for f in con.readline().split()[0:dim] ]) #Line 4: Box angles
        boxtemp=numpy.zeros((dim,dim),'d')
        boxtemp = length_angle_to_box(boxlengths,boxangles)
        con.readline() #Line 5: comment
        con.readline() #Line 6: comment
        num_types = int(con.readline().split()[0]) #Line 7: number of atom types
        num_each_type = con.readline().split() #line 8: number of each type of atom
        mass_of_type = con.readline().split() #line 9: mass of each type of atom
        num_atoms = 0
        for i in range(num_types):
            num_each_type[i] = int(num_each_type[i])
            mass_of_type[i] = float(mass_of_type[i])
            num_atoms += num_each_type[i]
        self.num_atoms = num_atoms
        self.box = boxtemp
        self.r = numpy.zeros((num_atoms,dim))
        self.free = numpy.zeros(num_atoms)
        #self.free = numpy.ones(num_atoms,dtype=bool)

        self.names = ['']*num_atoms
        self.mass = numpy.zeros(num_atoms)
        index = 0
        for i in range(num_types):
            name = con.readline().strip()
    #            if abs(1.0-mass_of_type[i]) < 1e-6 and name != "H":
    #                logger.warning("WARNING: Mass of %s set to 1.0", name)
            
            con.readline() #skip meaningless line
            for j in range(num_each_type[i]):
                vals = con.readline().split()
                for k in range(dim):
                    self.r[index][k] = float(vals[k])
                self.mass[index] = mass_of_type[i]
                self.names[index] = name
                if not int(vals[dim])==0:
                    #Think this should be maybe =1 ?
                    self.free[index]=0
                index += 1
    """

    def read_xyz(self, filename):
        '''
        reads a configuations from xyz file into the object
        Warning: destroys any previous content
        '''
        infile=open(filename,'r')
    	    	
        num_atoms=int( infile.readline().split()[0] )
   	
        self.r = numpy.zeros((num_atoms,3))
        self.free = numpy.ones(num_atoms)
        self.names = ['']*num_atoms
        self.mass = numpy.zeros(num_atoms)	
        self.num_atoms= num_atoms
        	
        tmpline=infile.readline().split()
        if(len(tmpline)==3):
            self.box[0][0]=tmpline[0]
            self.box[1][1]=tmpline[1]
            self.box[2][2]=tmpline[2]
    	for iatom in range(num_atoms):
    		tmpline=infile.readline().split()
    		self.names[iatom]=str(tmpline[0])
    		for k in range(3):
    			self.r[iatom][k]=float(tmpline[k+1])

    def write_xyz(self, filename=None):
        '''
        Load a con file
            filein: may be either a filename or a file-like object
        '''
        if filename==None:
            xyzfile=sys.stdout
        elif hasattr(filename, 'readline'):
            xyzfile=filename
        else:
            xyzfile = open(filename,'a')

        xyzfile.write("%d\n%f\t%f\t%f\n" % (self.num_atoms, self.box[0][0], self.box[1][1], self.box[2][2]) )
        for i in range(self.num_atoms):
            xyzfile.write("%s\t%f\t%f\t%f\n" % (self.names[i],self.r[i][0],self.r[i][1],self.r[i][2]) )

        if filename==None:
            pass
        elif hasattr(filename, 'readline'):
            pass
        else:
            xyzfile.close()

    def write_con(self,fileout,w='w'):
        '''
        Save a con file
            fileout: can be either a file name or a file-like object
            p:       information (in the form of an atoms object) to save
            w:       write/append flag
        '''
        if hasattr(fileout, 'write'):
            con = fileout
        else:
            con = open(fileout, w)
        print >> con, "Generated by eOn"
        print >> con
        dim = len(self.r[0])
        lengths, angles = box_to_length_angle(self.box)
        print >> con, " ".join(['%.4f' % s for s in lengths])
        print >> con, " ".join(['%.4f' % s for s in angles])
        print >> con
        print >> con
        
        type_count = [0]
        type_mass  = []
        type_name  = []
        
        type_name.append( self.names[0] )
        type_mass.append( self.mass[0] )
        
        for i in range(len(self)):
            name = self.names[i]
            #are we encountering a new type?
            if name != type_name[-1]:
                type_name.append(name)
                type_mass.append(self.mass[i])
                type_count.append(0)
            
            type_count[-1] +=1
        
        print >> con, len(type_name)
        print >> con, " ".join([str(i) for i in type_count])
        print >> con, " ".join(["%.4f"% i for i in type_mass])
        
        itype = 0
        iatom = 0
        for itype in range(len(type_name)):
            print >> con, type_name[itype]
            print >> con, "Coordinates of Component", itype+1
            for j in range(type_count[itype]):
                con.write("%.6f %.6f %.6f %d %d\n" %( self.r[iatom][0], self.r[iatom][1], self.r[iatom][2], int(not self.free[iatom]), iatom+1))
                iatom += 1

    """
    def write_conold(self,fileout, w = 'w'):
        '''
        Save a con file
            fileout: can be either a file name or a file-like object
            p:       information (in the form of an atoms object) to save
            w:       write/append flag
        '''
        if hasattr(fileout, 'write'):
            con = fileout
        else:
            con = open(fileout, w)
        print >> con, "Generated by eOn"
        print >> con
        dim = len(self.r[0])
        lengths, angles = box_to_length_angle(self.box)
        print >> con, " ".join(['%.6f' % s for s in lengths])
        print >> con, " ".join(['%.6f' % s for s in angles])
        print >> con
        print >> con
        atom_count = {}
        name_order = []
        for i in range(self.num_atoms):
            name = self.names[i]
            if name not in name_order:
                name_order.append(name)
            if name in atom_count:
                atom_count[name] += 1
            else:
                atom_count[name] = 1
        print >> con, len(name_order)
        print >> con, " ".join([str(atom_count[i]) for i in name_order])
        printmasses = []
        index = 0
        for i in range(len(name_order)):
            printmasses.append(self.mass[index])
            index += atom_count[name_order[i]]
        print >> con, " ".join(["%.6f"% i for i in printmasses])
        index = 0
        for i in range(len(name_order)):
            print >> con, name_order[i]
            print >> con, "Coordinates of Component", i+1
            for j in range(atom_count[name_order[i]]):
                con.write("%.6f %.6f %.6f %d %d\n" %( self.r[index][0], self.r[index][1], self.r[index][2], int(not self.free[index]), index+1))
                index += 1
    """
    def savecon(self,fileout, w = 'w'):
        '''
        Save a con file
            fileout: can be either a file name or a file-like object
            p:       information (in the form of an atoms object) to save
            w:       write/append flag
        '''
        if hasattr(fileout, 'write'):
            con = fileout
        else:
            con = open(fileout, w)
        print >> con, "Generated by eOn"
        print >> con
        dim = len(self.r[0])
        lengths, angles = box_to_length_angle(self.box)
        print >> con, " ".join(['%.4f' % s for s in lengths])
        print >> con, " ".join(['%.4f' % s for s in angles])
        print >> con
        print >> con
        
        type_count = [0]
        type_mass  = []
        type_name  = []
        
        type_name.append( self.names[0] )
        type_mass.append( self.mass[0] )
        
        for i in range(len(self)):
            name = self.names[i]
            #are we encountering a new type?
            if name != type_name[-1]:
                type_name.append(name)
                type_mass.append(self.mass[i])
                type_count.append(0)
            
            type_count[-1] +=1
    
        print >> con, len(type_name)
        print >> con, " ".join([str(i) for i in type_count])
        print >> con, " ".join(["%.4f"% i for i in type_mass])

        itype = 0
        iatom = 0
        for itype in range(len(type_name)):
            print >> con, type_name[itype]
            print >> con, "Coordinates of Component", itype+1
            for j in range(type_count[itype]):
                con.write("%.6f %.6f %.6f %d %d\n" %( self.r[iatom][0], self.r[iatom][1], self.r[iatom][2], int(not self.free[iatom]), iatom+1))
                iatom += 1



##########################################################################################

def pbc(r, box, ibox = None):
    """
    Applies periodic boundary conditions.
    Parameters:
        r:      the vector the boundary conditions are applied to
        box:    the box that defines the boundary conditions
        ibox:   the inverse of the box. This will be calcluated if not provided.
    """
    if ibox == None: 
        ibox = numpy.linalg.inv(box)
    vdir = numpy.dot(r, ibox)
    vdir = (vdir % 1.0 + 1.5) % 1.0 - 0.5
    return numpy.dot(vdir, box) 

def per_atom_norm(v, box, ibox = None):
    '''
    Returns a length N numpy array containing per atom distance
        v:      an Nx3 numpy array
        box:    box matrix that defines the boundary conditions
        ibox:   the inverse of the box. will be calculated if not provided
    '''
    diff = pbc(v, box, ibox)
    return numpy.array([numpy.linalg.norm(d) for d in diff])

def brute_neighbor_list(p, cutoff):
    nl = []
    ibox = numpy.linalg.inv(p.box)    
    for a in range(len(p)):
        nl.append([])
        for b in range(len(p)):
            if b != a:
                dist = numpy.linalg.norm(pbc(p.r[a] - p.r[b], p.box, ibox))        
                if dist < cutoff:
                    nl[a].append(b)
    return nl

def sweep_and_prune(p_in, cutoff, strict = True, bc = True):
    """ Returns a list of nearest neighbors within cutoff for each atom. 
        Parameters:
            p_in:   Atoms object
            cutoff: the radius within which two atoms are considered to intersect.
            strict: perform an actual distance check if True
            bc:     include neighbors across pbc's """
    #TODO: Get rid of 'cutoff' and use the covalent bond radii. (Rye can do)
        # Do we want to use covalent radii? I think the displace class wants to allow for user-defined cutoffs.
            # We should have both options available. -Rye
    #TODO: Make work for nonorthogonal boxes.
    p = p_in.copy()
    p.r = pbc(p.r, p.box)
    p.r -= numpy.array([min(p.r[:,0]), min(p.r[:,1]), min(p.r[:,2])]) 
    numatoms = len(p)
    coord_list = []
    for i in range(numatoms):
        coord_list.append([i, p.r[i]])
    for axis in range(3):
        sorted_axis = sorted(coord_list, key = lambda foo: foo[1][axis])
        intersect_axis = []
        for i in range(numatoms):
            intersect_axis.append([])
        for i in range(numatoms):
            done = False
            j = i + 1
            if not bc and j >= numatoms:
                done = True
            while not done:
                j = j % numatoms
                if j == i:
                    done = True
                dist = abs(sorted_axis[j][1][axis] - sorted_axis[i][1][axis])
                if p.box[axis][axis] - sorted_axis[i][1][axis] < cutoff:
                    dist = min(dist, (p.box[axis][axis] - sorted_axis[i][1][axis]) + sorted_axis[j][1][axis])
                if dist < cutoff:
                    intersect_axis[sorted_axis[i][0]].append(sorted_axis[j][0]) 
                    intersect_axis[sorted_axis[j][0]].append(sorted_axis[i][0]) 
                    j += 1
                    if not bc and j >= numatoms:
                        done = True
                else:
                    done = True
        if axis == 0:
            intersect = []
            for i in range(numatoms):
                intersect.append([])
                intersect[i] = intersect_axis[i]
        else:
            for i in range(numatoms):
                intersect[i] = list(set(intersect[i]).intersection(intersect_axis[i]))
    if strict:
        ibox = numpy.linalg.inv(p.box)
        for i in range(numatoms):
            l = intersect[i][:]
            for j in l:
                dist = numpy.linalg.norm(pbc(p.r[i] - p.r[j], p.box, ibox))
                if dist > cutoff:
                    intersect[i].remove(j)
                    intersect[j].remove(i)
    return intersect

def neighbor_list(p, cutoff, brute=False):
    if brute:
        nl = brute_neighbor_list(p, cutoff)
    else:
        nl = sweep_and_prune(p, cutoff)
    return nl

def box_to_length_angle(box):
    lengths = numpy.zeros(3)
    lengths[0] = numpy.linalg.norm(box[0,:])
    lengths[1] = numpy.linalg.norm(box[1,:])
    lengths[2] = numpy.linalg.norm(box[2,:])
    angles = numpy.zeros(3)
    angles[0] = numpy.arccos(numpy.dot(box[0,:]/lengths[0],box[1,:]/lengths[1]))
    angles[1] = numpy.arccos(numpy.dot(box[0,:]/lengths[0],box[2,:]/lengths[2]))
    angles[2] = numpy.arccos(numpy.dot(box[1,:]/lengths[1],box[2,:]/lengths[2]))
    angles *= 180.0/numpy.pi
    return lengths, angles


def match(a,b,eps_r):
    if len(a)!=len(b):
        return False
    
    return max(per_atom_norm(a.r-b.r, a.box))<eps_r




elements = {}
numElements = 119
elements[  0] = elements[ 'Xx'] = {'symbol':  'Xx', 'name':       'unknown', 'mass':   1.00000000, 'radius':  1.0000, 'color': [1.000, 0.078, 0.576], 'number': 0}
elements[  1] = elements[  'H'] = {'symbol':   'H', 'name':      'hydrogen', 'mass':   1.00794000, 'radius':  0.3100, 'color': [1.000, 1.000, 1.000], 'number': 1}
elements[  2] = elements[ 'He'] = {'symbol':  'He', 'name':        'helium', 'mass':   4.00260200, 'radius':  0.2800, 'color': [0.851, 1.000, 1.000], 'number': 2}
elements[  3] = elements[ 'Li'] = {'symbol':  'Li', 'name':       'lithium', 'mass':   6.94100000, 'radius':  1.2800, 'color': [0.800, 0.502, 1.000], 'number': 3}
elements[  4] = elements[ 'Be'] = {'symbol':  'Be', 'name':     'beryllium', 'mass':   9.01218200, 'radius':  0.9600, 'color': [0.761, 1.000, 0.000], 'number': 4}
elements[  5] = elements[  'B'] = {'symbol':   'B', 'name':         'boron', 'mass':  10.81100000, 'radius':  0.8400, 'color': [1.000, 0.710, 0.710], 'number': 5}
elements[  6] = elements[  'C'] = {'symbol':   'C', 'name':        'carbon', 'mass':  12.01070000, 'radius':  0.7300, 'color': [0.565, 0.565, 0.565], 'number': 6}
elements[  7] = elements[  'N'] = {'symbol':   'N', 'name':      'nitrogen', 'mass':  14.00670000, 'radius':  0.7100, 'color': [0.188, 0.314, 0.973], 'number': 7}
elements[  8] = elements[  'O'] = {'symbol':   'O', 'name':        'oxygen', 'mass':  15.99940000, 'radius':  0.6600, 'color': [1.000, 0.051, 0.051], 'number': 8}
elements[  9] = elements[  'F'] = {'symbol':   'F', 'name':      'fluorine', 'mass':  18.99840320, 'radius':  0.5700, 'color': [0.565, 0.878, 0.314], 'number': 9}
elements[ 10] = elements[ 'Ne'] = {'symbol':  'Ne', 'name':          'neon', 'mass':  20.17970000, 'radius':  0.5800, 'color': [0.702, 0.890, 0.961], 'number': 10}
elements[ 11] = elements[ 'Na'] = {'symbol':  'Na', 'name':        'sodium', 'mass':  22.98976928, 'radius':  1.6600, 'color': [0.671, 0.361, 0.949], 'number': 11}
elements[ 12] = elements[ 'Mg'] = {'symbol':  'Mg', 'name':     'magnesium', 'mass':  24.30500000, 'radius':  1.4100, 'color': [0.541, 1.000, 0.000], 'number': 12}
elements[ 13] = elements[ 'Al'] = {'symbol':  'Al', 'name':      'aluminum', 'mass':  26.98153860, 'radius':  1.2100, 'color': [0.749, 0.651, 0.651], 'number': 13}
elements[ 14] = elements[ 'Si'] = {'symbol':  'Si', 'name':       'silicon', 'mass':  28.08550000, 'radius':  1.1100, 'color': [0.941, 0.784, 0.627], 'number': 14}
elements[ 15] = elements[  'P'] = {'symbol':   'P', 'name':    'phosphorus', 'mass':  30.97376200, 'radius':  1.0700, 'color': [1.000, 0.502, 0.000], 'number': 15}
elements[ 16] = elements[  'S'] = {'symbol':   'S', 'name':        'sulfur', 'mass':  32.06500000, 'radius':  1.0500, 'color': [1.000, 1.000, 0.188], 'number': 16}
elements[ 17] = elements[ 'Cl'] = {'symbol':  'Cl', 'name':      'chlorine', 'mass':  35.45300000, 'radius':  1.0200, 'color': [0.122, 0.941, 0.122], 'number': 17}
elements[ 18] = elements[ 'Ar'] = {'symbol':  'Ar', 'name':         'argon', 'mass':  39.94800000, 'radius':  1.0600, 'color': [0.502, 0.820, 0.890], 'number': 18}
elements[ 19] = elements[  'K'] = {'symbol':   'K', 'name':     'potassium', 'mass':  39.09830000, 'radius':  2.0300, 'color': [0.561, 0.251, 0.831], 'number': 19}
elements[ 20] = elements[ 'Ca'] = {'symbol':  'Ca', 'name':       'calcium', 'mass':  40.07800000, 'radius':  1.7600, 'color': [0.239, 1.000, 0.000], 'number': 20}
elements[ 21] = elements[ 'Sc'] = {'symbol':  'Sc', 'name':      'scandium', 'mass':  44.95591200, 'radius':  1.7000, 'color': [0.902, 0.902, 0.902], 'number': 21}
elements[ 22] = elements[ 'Ti'] = {'symbol':  'Ti', 'name':      'titanium', 'mass':  47.86700000, 'radius':  1.6000, 'color': [0.749, 0.761, 0.780], 'number': 22}
elements[ 23] = elements[  'V'] = {'symbol':   'V', 'name':      'vanadium', 'mass':  50.94150000, 'radius':  1.5300, 'color': [0.651, 0.651, 0.671], 'number': 23}
elements[ 24] = elements[ 'Cr'] = {'symbol':  'Cr', 'name':      'chromium', 'mass':  51.99610000, 'radius':  1.3900, 'color': [0.541, 0.600, 0.780], 'number': 24}
elements[ 25] = elements[ 'Mn'] = {'symbol':  'Mn', 'name':     'manganese', 'mass':  54.93804500, 'radius':  1.3900, 'color': [0.611, 0.478, 0.780], 'number': 25}
elements[ 26] = elements[ 'Fe'] = {'symbol':  'Fe', 'name':          'iron', 'mass':  55.84500000, 'radius':  1.3200, 'color': [0.878, 0.400, 0.200], 'number': 26}
elements[ 27] = elements[ 'Co'] = {'symbol':  'Co', 'name':        'cobalt', 'mass':  58.69340000, 'radius':  1.2600, 'color': [0.941, 0.565, 0.627], 'number': 27}
elements[ 28] = elements[ 'Ni'] = {'symbol':  'Ni', 'name':        'nickel', 'mass':  58.93319500, 'radius':  1.2400, 'color': [0.314, 0.816, 0.314], 'number': 28}
elements[ 29] = elements[ 'Cu'] = {'symbol':  'Cu', 'name':        'copper', 'mass':  63.54600000, 'radius':  1.3200, 'color': [0.784, 0.502, 0.200], 'number': 29}
elements[ 30] = elements[ 'Zn'] = {'symbol':  'Zn', 'name':          'zinc', 'mass':  65.38000000, 'radius':  1.2200, 'color': [0.490, 0.502, 0.690], 'number': 30}
elements[ 31] = elements[ 'Ga'] = {'symbol':  'Ga', 'name':       'gallium', 'mass':  69.72300000, 'radius':  1.2200, 'color': [0.761, 0.561, 0.561], 'number': 31}
elements[ 32] = elements[ 'Ge'] = {'symbol':  'Ge', 'name':     'germanium', 'mass':  72.64000000, 'radius':  1.2000, 'color': [0.400, 0.561, 0.561], 'number': 32}
elements[ 33] = elements[ 'As'] = {'symbol':  'As', 'name':       'arsenic', 'mass':  74.92160000, 'radius':  1.1900, 'color': [0.741, 0.502, 0.890], 'number': 33}
elements[ 34] = elements[ 'Se'] = {'symbol':  'Se', 'name':      'selenium', 'mass':  78.96000000, 'radius':  1.2000, 'color': [1.000, 0.631, 0.000], 'number': 34}
elements[ 35] = elements[ 'Br'] = {'symbol':  'Br', 'name':       'bromine', 'mass':  79.90400000, 'radius':  1.2000, 'color': [0.651, 0.161, 0.161], 'number': 35}
elements[ 36] = elements[ 'Kr'] = {'symbol':  'Kr', 'name':       'krypton', 'mass':  83.79800000, 'radius':  1.1600, 'color': [0.361, 0.722, 0.820], 'number': 36}
elements[ 37] = elements[ 'Rb'] = {'symbol':  'Rb', 'name':      'rubidium', 'mass':  85.46780000, 'radius':  2.2000, 'color': [0.439, 0.180, 0.690], 'number': 37}
elements[ 38] = elements[ 'Sr'] = {'symbol':  'Sr', 'name':     'strontium', 'mass':  87.62000000, 'radius':  1.9500, 'color': [0.000, 1.000, 0.000], 'number': 38}
elements[ 39] = elements[  'Y'] = {'symbol':   'Y', 'name':       'yttrium', 'mass':  88.90585000, 'radius':  1.9000, 'color': [0.580, 1.000, 1.000], 'number': 39}
elements[ 40] = elements[ 'Zr'] = {'symbol':  'Zr', 'name':     'zirconium', 'mass':  91.22400000, 'radius':  1.7500, 'color': [0.580, 0.878, 0.878], 'number': 40}
elements[ 41] = elements[ 'Nb'] = {'symbol':  'Nb', 'name':       'niobium', 'mass':  92.90638000, 'radius':  1.6400, 'color': [0.451, 0.761, 0.788], 'number': 41}
elements[ 42] = elements[ 'Mo'] = {'symbol':  'Mo', 'name':    'molybdenum', 'mass':  95.96000000, 'radius':  1.5400, 'color': [0.329, 0.710, 0.710], 'number': 42}
elements[ 43] = elements[ 'Tc'] = {'symbol':  'Tc', 'name':    'technetium', 'mass':  98.00000000, 'radius':  1.4700, 'color': [0.231, 0.620, 0.620], 'number': 43}
elements[ 44] = elements[ 'Ru'] = {'symbol':  'Ru', 'name':     'ruthenium', 'mass': 101.07000000, 'radius':  1.4600, 'color': [0.141, 0.561, 0.561], 'number': 44}
elements[ 45] = elements[ 'Rh'] = {'symbol':  'Rh', 'name':       'rhodium', 'mass': 102.90550000, 'radius':  1.4200, 'color': [0.039, 0.490, 0.549], 'number': 45}
elements[ 46] = elements[ 'Pd'] = {'symbol':  'Pd', 'name':     'palladium', 'mass': 106.42000000, 'radius':  1.3900, 'color': [0.000, 0.412, 0.522], 'number': 46}
elements[ 47] = elements[ 'Ag'] = {'symbol':  'Ag', 'name':        'silver', 'mass': 107.86820000, 'radius':  1.4500, 'color': [0.753, 0.753, 0.753], 'number': 47}
elements[ 48] = elements[ 'Cd'] = {'symbol':  'Cd', 'name':       'cadmium', 'mass': 112.41100000, 'radius':  1.4400, 'color': [1.000, 0.851, 0.561], 'number': 48}
elements[ 49] = elements[ 'In'] = {'symbol':  'In', 'name':        'indium', 'mass': 114.81800000, 'radius':  1.4200, 'color': [0.651, 0.459, 0.451], 'number': 49}
elements[ 50] = elements[ 'Sn'] = {'symbol':  'Sn', 'name':           'tin', 'mass': 118.71000000, 'radius':  1.3900, 'color': [0.400, 0.502, 0.502], 'number': 50}
elements[ 51] = elements[ 'Sb'] = {'symbol':  'Sb', 'name':      'antimony', 'mass': 121.76000000, 'radius':  1.3900, 'color': [0.620, 0.388, 0.710], 'number': 51}
elements[ 52] = elements[ 'Te'] = {'symbol':  'Te', 'name':     'tellurium', 'mass': 127.60000000, 'radius':  1.3800, 'color': [0.831, 0.478, 0.000], 'number': 52}
elements[ 53] = elements[  'I'] = {'symbol':   'I', 'name':        'iodine', 'mass': 126.90470000, 'radius':  1.3900, 'color': [0.580, 0.000, 0.580], 'number': 53}
elements[ 54] = elements[ 'Xe'] = {'symbol':  'Xe', 'name':         'xenon', 'mass': 131.29300000, 'radius':  1.4000, 'color': [0.259, 0.620, 0.690], 'number': 54}
elements[ 55] = elements[ 'Cs'] = {'symbol':  'Cs', 'name':        'cesium', 'mass': 132.90545190, 'radius':  2.4400, 'color': [0.341, 0.090, 0.561], 'number': 55}
elements[ 56] = elements[ 'Ba'] = {'symbol':  'Ba', 'name':        'barium', 'mass': 137.32700000, 'radius':  2.1500, 'color': [0.000, 0.788, 0.000], 'number': 56}
elements[ 57] = elements[ 'La'] = {'symbol':  'La', 'name':     'lanthanum', 'mass': 138.90547000, 'radius':  2.0700, 'color': [0.439, 0.831, 1.000], 'number': 57}
elements[ 58] = elements[ 'Ce'] = {'symbol':  'Ce', 'name':        'cerium', 'mass': 140.11600000, 'radius':  2.0400, 'color': [1.000, 1.000, 0.780], 'number': 58}
elements[ 59] = elements[ 'Pr'] = {'symbol':  'Pr', 'name':  'praseodymium', 'mass': 140.90765000, 'radius':  2.0300, 'color': [0.851, 1.000, 0.780], 'number': 59}
elements[ 60] = elements[ 'Nd'] = {'symbol':  'Nd', 'name':     'neodymium', 'mass': 144.24200000, 'radius':  2.0100, 'color': [0.780, 1.000, 0.780], 'number': 60}
elements[ 61] = elements[ 'Pm'] = {'symbol':  'Pm', 'name':    'promethium', 'mass': 145.00000000, 'radius':  1.9900, 'color': [0.639, 1.000, 0.780], 'number': 61}
elements[ 62] = elements[ 'Sm'] = {'symbol':  'Sm', 'name':      'samarium', 'mass': 150.36000000, 'radius':  1.9800, 'color': [0.561, 1.000, 0.780], 'number': 62}
elements[ 63] = elements[ 'Eu'] = {'symbol':  'Eu', 'name':      'europium', 'mass': 151.96400000, 'radius':  1.9800, 'color': [0.380, 1.000, 0.780], 'number': 63}
elements[ 64] = elements[ 'Gd'] = {'symbol':  'Gd', 'name':    'gadolinium', 'mass': 157.25000000, 'radius':  1.9600, 'color': [0.271, 1.000, 0.780], 'number': 64}
elements[ 65] = elements[ 'Tb'] = {'symbol':  'Tb', 'name':       'terbium', 'mass': 158.92535000, 'radius':  1.9400, 'color': [0.189, 1.000, 0.780], 'number': 65}
elements[ 66] = elements[ 'Dy'] = {'symbol':  'Dy', 'name':    'dysprosium', 'mass': 162.50000000, 'radius':  1.9200, 'color': [0.122, 1.000, 0.780], 'number': 66}
elements[ 67] = elements[ 'Ho'] = {'symbol':  'Ho', 'name':       'holmium', 'mass': 164.93032000, 'radius':  1.9200, 'color': [0.000, 1.000, 0.612], 'number': 67}
elements[ 68] = elements[ 'Er'] = {'symbol':  'Er', 'name':        'erbium', 'mass': 167.25900000, 'radius':  1.8900, 'color': [0.000, 0.902, 0.459], 'number': 68}
elements[ 69] = elements[ 'Tm'] = {'symbol':  'Tm', 'name':       'thulium', 'mass': 168.93421000, 'radius':  1.9000, 'color': [0.000, 0.831, 0.322], 'number': 69}
elements[ 70] = elements[ 'Yb'] = {'symbol':  'Yb', 'name':     'ytterbium', 'mass': 173.05400000, 'radius':  1.8700, 'color': [0.000, 0.749, 0.220], 'number': 70}
elements[ 71] = elements[ 'Lu'] = {'symbol':  'Lu', 'name':      'lutetium', 'mass': 174.96680000, 'radius':  1.8700, 'color': [0.000, 0.671, 0.141], 'number': 71}
elements[ 72] = elements[ 'Hf'] = {'symbol':  'Hf', 'name':       'hafnium', 'mass': 178.49000000, 'radius':  1.7500, 'color': [0.302, 0.761, 1.000], 'number': 72}
elements[ 73] = elements[ 'Ta'] = {'symbol':  'Ta', 'name':      'tantalum', 'mass': 180.94788000, 'radius':  1.7000, 'color': [0.302, 0.651, 1.000], 'number': 73}
elements[ 74] = elements[  'W'] = {'symbol':   'W', 'name':      'tungsten', 'mass': 183.84000000, 'radius':  1.6200, 'color': [0.129, 0.580, 0.839], 'number': 74}
elements[ 75] = elements[ 'Re'] = {'symbol':  'Re', 'name':       'rhenium', 'mass': 186.20700000, 'radius':  1.5100, 'color': [0.149, 0.490, 0.671], 'number': 75}
elements[ 76] = elements[ 'Os'] = {'symbol':  'Os', 'name':        'osmium', 'mass': 190.23000000, 'radius':  1.4400, 'color': [0.149, 0.400, 0.588], 'number': 76}
elements[ 77] = elements[ 'Ir'] = {'symbol':  'Ir', 'name':       'iridium', 'mass': 192.21700000, 'radius':  1.4100, 'color': [0.090, 0.329, 0.529], 'number': 77}
elements[ 78] = elements[ 'Pt'] = {'symbol':  'Pt', 'name':      'platinum', 'mass': 195.08400000, 'radius':  1.3600, 'color': [0.816, 0.816, 0.878], 'number': 78}
elements[ 79] = elements[ 'Au'] = {'symbol':  'Au', 'name':          'gold', 'mass': 196.96656900, 'radius':  1.3600, 'color': [1.000, 0.820, 0.137], 'number': 79}
elements[ 80] = elements[ 'Hg'] = {'symbol':  'Hg', 'name':       'mercury', 'mass': 200.59000000, 'radius':  1.3200, 'color': [0.722, 0.722, 0.816], 'number': 80}
elements[ 81] = elements[ 'Tl'] = {'symbol':  'Tl', 'name':      'thallium', 'mass': 204.38330000, 'radius':  1.4500, 'color': [0.651, 0.329, 0.302], 'number': 81}
elements[ 82] = elements[ 'Pb'] = {'symbol':  'Pb', 'name':          'lead', 'mass': 207.20000000, 'radius':  1.4600, 'color': [0.341, 0.349, 0.380], 'number': 82}
elements[ 83] = elements[ 'Bi'] = {'symbol':  'Bi', 'name':       'bismuth', 'mass': 208.98040000, 'radius':  1.4800, 'color': [0.620, 0.310, 0.710], 'number': 83}
elements[ 84] = elements[ 'Po'] = {'symbol':  'Po', 'name':      'polonium', 'mass': 210.00000000, 'radius':  1.4000, 'color': [0.671, 0.361, 0.000], 'number': 84}
elements[ 85] = elements[ 'At'] = {'symbol':  'At', 'name':      'astatine', 'mass': 210.00000000, 'radius':  1.5000, 'color': [0.459, 0.310, 0.271], 'number': 85}
elements[ 86] = elements[ 'Rn'] = {'symbol':  'Rn', 'name':         'radon', 'mass': 220.00000000, 'radius':  1.5000, 'color': [0.259, 0.510, 0.588], 'number': 86}
elements[ 87] = elements[ 'Fr'] = {'symbol':  'Fr', 'name':      'francium', 'mass': 223.00000000, 'radius':  2.6000, 'color': [0.259, 0.000, 0.400], 'number': 87}
elements[ 88] = elements[ 'Ra'] = {'symbol':  'Ra', 'name':        'radium', 'mass': 226.00000000, 'radius':  2.2100, 'color': [0.000, 0.490, 0.000], 'number': 88}
elements[ 89] = elements[ 'Ac'] = {'symbol':  'Ac', 'name':      'actinium', 'mass': 227.00000000, 'radius':  2.1500, 'color': [0.439, 0.671, 0.980], 'number': 89}
elements[ 90] = elements[ 'Th'] = {'symbol':  'Th', 'name':       'thorium', 'mass': 231.03588000, 'radius':  2.0600, 'color': [0.000, 0.729, 1.000], 'number': 90}
elements[ 91] = elements[ 'Pa'] = {'symbol':  'Pa', 'name':  'protactinium', 'mass': 232.03806000, 'radius':  2.0000, 'color': [0.000, 0.631, 1.000], 'number': 91}
elements[ 92] = elements[  'U'] = {'symbol':   'U', 'name':       'uranium', 'mass': 237.00000000, 'radius':  1.9600, 'color': [0.000, 0.561, 1.000], 'number': 92}
elements[ 93] = elements[ 'Np'] = {'symbol':  'Np', 'name':     'neptunium', 'mass': 238.02891000, 'radius':  1.9000, 'color': [0.000, 0.502, 1.000], 'number': 93}
elements[ 94] = elements[ 'Pu'] = {'symbol':  'Pu', 'name':     'plutonium', 'mass': 243.00000000, 'radius':  1.8700, 'color': [0.000, 0.420, 1.000], 'number': 94}
elements[ 95] = elements[ 'Am'] = {'symbol':  'Am', 'name':     'americium', 'mass': 244.00000000, 'radius':  1.8000, 'color': [0.329, 0.361, 0.949], 'number': 95}
elements[ 96] = elements[ 'Cm'] = {'symbol':  'Cm', 'name':        'curium', 'mass': 247.00000000, 'radius':  1.6900, 'color': [0.471, 0.361, 0.890], 'number': 96}
elements[ 97] = elements[ 'Bk'] = {'symbol':  'Bk', 'name':     'berkelium', 'mass': 247.00000000, 'radius':  1.6600, 'color': [0.541, 0.310, 0.890], 'number': 97}
elements[ 98] = elements[ 'Cf'] = {'symbol':  'Cf', 'name':   'californium', 'mass': 251.00000000, 'radius':  1.6800, 'color': [0.631, 0.212, 0.831], 'number': 98}
elements[ 99] = elements[ 'Es'] = {'symbol':  'Es', 'name':   'einsteinium', 'mass': 252.00000000, 'radius':  1.6500, 'color': [0.702, 0.122, 0.831], 'number': 99}
elements[100] = elements[ 'Fm'] = {'symbol':  'Fm', 'name':       'fermium', 'mass': 257.00000000, 'radius':  1.6700, 'color': [0.702, 0.122, 0.729], 'number': 100}
elements[101] = elements[ 'Md'] = {'symbol':  'Md', 'name':   'mendelevium', 'mass': 258.00000000, 'radius':  1.7300, 'color': [0.702, 0.051, 0.651], 'number': 101}
elements[102] = elements[ 'No'] = {'symbol':  'No', 'name':      'nobelium', 'mass': 259.00000000, 'radius':  1.7600, 'color': [0.741, 0.051, 0.529], 'number': 102}
elements[103] = elements[ 'Lr'] = {'symbol':  'Lr', 'name':    'lawrencium', 'mass': 262.00000000, 'radius':  1.6100, 'color': [0.780, 0.000, 0.400], 'number': 103}
elements[104] = elements[ 'Rf'] = {'symbol':  'Rf', 'name': 'rutherfordium', 'mass': 261.00000000, 'radius':  1.5700, 'color': [0.800, 0.000, 0.349], 'number': 104}
elements[105] = elements[ 'Db'] = {'symbol':  'Db', 'name':       'dubnium', 'mass': 262.00000000, 'radius':  1.4900, 'color': [0.820, 0.000, 0.310], 'number': 105}
elements[106] = elements[ 'Sg'] = {'symbol':  'Sg', 'name':    'seaborgium', 'mass': 266.00000000, 'radius':  1.4300, 'color': [0.851, 0.000, 0.271], 'number': 106}
elements[107] = elements[ 'Bh'] = {'symbol':  'Bh', 'name':       'bohrium', 'mass': 264.00000000, 'radius':  1.4100, 'color': [0.878, 0.000, 0.220], 'number': 107}
elements[108] = elements[ 'Hs'] = {'symbol':  'Hs', 'name':       'hassium', 'mass': 277.00000000, 'radius':  1.3400, 'color': [0.902, 0.000, 0.180], 'number': 108}
elements[109] = elements[ 'Mt'] = {'symbol':  'Mt', 'name':    'meitnerium', 'mass': 268.00000000, 'radius':  1.2900, 'color': [0.922, 0.000, 0.149], 'number': 109}
elements[110] = elements[ 'Ds'] = {'symbol':  'Ds', 'name':            'Ds', 'mass': 271.00000000, 'radius':  1.2800, 'color': [0.922, 0.000, 0.149], 'number': 110}
elements[111] = elements['Uuu'] = {'symbol': 'Uuu', 'name':           'Uuu', 'mass': 272.00000000, 'radius':  1.2100, 'color': [0.922, 0.000, 0.149], 'number': 111}
elements[112] = elements['Uub'] = {'symbol': 'Uub', 'name':           'Uub', 'mass': 285.00000000, 'radius':  1.2200, 'color': [0.922, 0.000, 0.149], 'number': 112}
elements[113] = elements['Uut'] = {'symbol': 'Uut', 'name':           'Uut', 'mass': 284.00000000, 'radius':  1.3600, 'color': [0.922, 0.000, 0.149], 'number': 113}
elements[114] = elements['Uuq'] = {'symbol': 'Uuq', 'name':           'Uuq', 'mass': 289.00000000, 'radius':  1.4300, 'color': [0.922, 0.000, 0.149], 'number': 114}
elements[115] = elements['Uup'] = {'symbol': 'Uup', 'name':           'Uup', 'mass': 288.00000000, 'radius':  1.6200, 'color': [0.922, 0.000, 0.149], 'number': 115}
elements[116] = elements['Uuh'] = {'symbol': 'Uuh', 'name':           'Uuh', 'mass': 292.00000000, 'radius':  1.7500, 'color': [0.922, 0.000, 0.149], 'number': 116}
elements[117] = elements['Uus'] = {'symbol': 'Uus', 'name':           'Uus', 'mass': 294.00000000, 'radius':  1.6500, 'color': [0.922, 0.000, 0.149], 'number': 117}
elements[118] = elements['Uuo'] = {'symbol': 'Uuo', 'name':           'Uuo', 'mass': 296.00000000, 'radius':  1.5700, 'color': [0.922, 0.000, 0.149], 'number': 118}


