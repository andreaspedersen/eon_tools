"""
    Author: Kjartan Myrdal and others
    Date: June 2012
    License: Undefined, if you have this you're probably allowed to use it.
"""

import numpy as np
try:
    from ase.data import atomic_numbers
    import ase
    isAse = True
except ImportError:
    print 'Ase module not found.'
    isAse = False
import eon_tools as et

_floatErr = 1e-6


#from ase import Atom
#from ase import Atoms
from scipy import zeros, array, floor
import re
import atom



def elements_in_line(line):
    element = ''
    elements = list()
    lastchar = ' '
    index = 0
    for char in line:
        if char != ' ' and char != '\n' and char !='\t':
            if lastchar == ' ' or lastchar =='\t':
                element=char
            else:
                element=element+char
        elif lastchar != ' ' and lastchar !='\t':
            elements.append(element)
        lastchar = char
    return elements

def is_number(n):
    if re.match("^[-+]?[0-9]+$", n):
        return True
    return False

def _unkownElement(element):
    print "Element %s is unknown." %(element,)
    return 1


def con_2_eon_atoms(path):
    """
        Loads con file
        Takes in a path to a .con file and returns
        eon_atoms object.
        """
    at = et.Atoms(path)
    #    at = et.Atoms(0)
    #    at.read_con(path)
    return at

def con_2_array(fileName = None, boundary = [1,1,1]):
    if fileName == None:
        fileName = raw_input("Enter file name (.con is added): ")
        fileName = fileName+".con"
    
    f = open(fileName,'r')
	
    ## seed to random
    line = f.readline()
    ## time stuff
    line = f.readline()
    ## box size
    box = f.readline().split()
    box[0] = float(box[0])
    box[1] = float(box[1])
    box[2] = float(box[2])
    
    line = f.readline()
    angle_basis_vec = elements_in_line(line)
    if int(float(angle_basis_vec[0])) != 90 or int(float(angle_basis_vec[1])) != 90 or int(float(angle_basis_vec[2])) != 90:
        print "Only orthogonal lattice vectores are accepted!"
        atoms = None
	
    else:
        ## pressure stuff
        line = f.readline()
        ## stuff
        line = f.readline()
		## number of different elements
        line = f.readline()
        nrTypes = int(elements_in_line(line)[0])
        ## number of atoms of the different kinds
        line = f.readline()
        nrAtomsPerType = elements_in_line(line)
        nrAtoms=0
        for nr in nrAtomsPerType:
            nrAtoms=nrAtoms+int(nr)
        ## the different types mass
        line = f.readline()
        
        atomPos=zeros([nrAtoms,3],'f')
        for i in range(0,nrTypes):
            ## type of atom
            line = f.readline()
            index=0
            
            ## text telling the components number
            line = f.readline()
            for j in range(0,int(nrAtomsPerType[i])):
                ## info about a certain atom
                line = f.readline()
                atomInfo = elements_in_line(line)
                
                atomPos[index,0]=float(atomInfo[0])
                atomPos[index,1]=float(atomInfo[1])
                atomPos[index,2]=float(atomInfo[2])
                
                index=index+1
    f.close()
    return (box, atomPos)


if isAse:
## ASE HELPER STUFF START
    
    def diff(start,end):
        unitCell = end.get_cell()
        unitCellVectors = array([unitCell[0,0], unitCell[1,1], unitCell[2,2]])
        diffS = end.get_scaled_positions()-start.get_scaled_positions()
        diffS = diffS-floor(diffS+0.5)
        diffC = diffS*unitCellVectors    
        return (diffC,diffS) 
        
    def correct_pos(atoms):
        unitCell = atoms.get_cell()
        unitCellVectors = array([unitCell[0,0], unitCell[1,1], unitCell[2,2]])
        posAll = atoms.get_positions()
        ## correct x
        posx = unitCellVectors[0]*(posAll[:,0] < 0)-unitCellVectors[0]*(unitCellVectors[0] < posAll[:,0])
        ## correct y
        posy = unitCellVectors[1]*(posAll[:,1] < 0)-unitCellVectors[1]*(unitCellVectors[1] < posAll[:,1])
        ## correct z
        posz = unitCellVectors[2]*(posAll[:,2] < 0)-unitCellVectors[2]*(unitCellVectors[2] < posAll[:,2])
        
        for i in range(0,len(posx)):
            posAll[i] = posAll[i]+array([posx[i],posy[i],posz[i]])
        atoms.set_positions(posAll)
        return
        
    def move_2_origo(atoms):
        pos = atoms.GetCartesianPositions()
        x = min(pos[:,0])
        y = min(pos[:,1])
        z = min(pos[:,2])
        atoms.SetCartesianPositions(pos-[x,y,z]+[0.1,0.1,0.1])
        return
        
    def displace_atoms(atoms, displacement):
        pos = atoms.get_positions()
        atoms.set_positions(pos+array(displacement))
        correctPos(atoms)
        return

    def sort_types(atoms):
        atomic_nrs = atoms.get_atomic_numbers()
        types = []
        for atomic_nr in atomic_nrs:
            if atomic_nr not in types:
                types.append(atomic_nr)

        types.sort()

        for type in types:
            #sort atoms according to type:
            for i in range(len(atoms)-1, -1, -1):
                if atomic_nrs[i] == type:
                    atoms.append(atoms.pop(i))
            atomic_nrs = atoms.get_atomic_numbers()
        return

    def ase_change_atom_nr(atoms, nr_old, nr_new):
        nrs = atoms.get_atomic_numbers()
        
        for i in range(0, len(nrs)):
            if nrs[i] == nr_old:
                nrs[i] = nr_new
        
        atoms.set_atomic_numbers(nrs)
        return

    def ase_translate_and_wrap(atoms,disp=[10,10,10]):
        atoms.translate(disp)
        atoms.set_scaled_positions(atoms.get_scaled_positions())
        return

## ASE HELPER STUFF END
##-----------------------
## ASE FROM AND TO SOMETHING START

    def con_2_ase(path, boundary = [1,1,1]):
        
        if not path.endswith('.con'):
            path += '.con'
        f = open(path, 'r')
        
        ## seed to random
        f.readline()
        
        ## time stuff
        f.readline()
        
        basis_vec = map(float, f.readline().split())
        
        angle_basis_vec = map(lambda x: x - 90 < _floatErr, map(float, f.readline().split()))
        if not np.array(angle_basis_vec).sum():
            print "Only orthogonal lattice vectores are accepted!"
            return
        
        f.readline()
        f.readline()
        nrComponents = int(f.readline().strip())
        nrAtomsPrType = np.array(map(int, f.readline().split()))
        
        #Masses, skiped in original parser. Same here
        f.readline()
        
        #All data is parced into these arrays before Atoms object is created.
        pos = np.zeros((nrAtomsPrType.sum(), 3))
        atomNumbers = np.zeros(nrAtomsPrType.sum())
        filt = np.zeros(nrAtomsPrType.sum(), dtype='int32')
        
        #
        for comp in xrange(nrComponents):
            element = f.readline().strip()
            if element.isdigit():
                element = int(element)
            else:
                element = atomic_numbers.get(element) or self._unkownElement(element)
            
            f.readline()
            
            for i in xrange(nrAtomsPrType[0:comp].sum(), nrAtomsPrType[0:comp+1].sum()):
                data = f.readline().split()
                pos[i] = map(float, data[:3])
                atomNumbers[i] = element
                filt[i] = int(data[3])
        
        f.close()
        return (ase.Atoms(numbers=atomNumbers, positions=pos, cell=basis_vec, pbc=boundary), abs(filt-1))


    # convert atoms object
    def eon_2_ase(e_atoms, boundary = [1,1,1]):
        return ase.Atoms(numbers=e_atoms.get_atomic_numbers(), positions=e_atoms.r.copy(), cell=e_atoms.box.copy(), pbc=boundary)

    def lammps_xyz_2_ase(fileName, boundary = [1,1,1]):
        f = open(fileName,'r')
        line = f.readline()
        line = f.readline()
        line = f.readline()
        line = f.readline()
        data = elements_in_line(line)
        atomNrs = int(data[0])
        
        line = f.readline()
        line = f.readline()
        data = elements_in_line(line)
        xMin = float(data[0])
        xMax = float(data[1])
        u1 = xMax-xMin
        
        line = f.readline()
        data = elements_in_line(line)
        yMin = float(data[0])
        yMax = float(data[1])
        u2 = yMax-yMin
        
        line = f.readline()
        data = elements_in_line(line)
        zMin = float(data[0])
        zMax = float(data[1])
        u3 = zMax-zMin
        
        atoms_read = {}
        
        line = f.readline()
        for i in range(atomNrs):
            line = f.readline()
            data = elements_in_line(line)
            atoms_read[int(data[0])-1] = [int(data[1]), [float(data[2]), float(data[3]), float(data[4])]]
        
        atomlist = []
        
        keys = atoms_read.keys()
        keys.sort()
        for key in keys:
            atomlist.append(ase.Atom( symbol=atoms_read[key][0], position=atoms_read[key][1] ))
        
        atoms = ase.Atoms(atomlist,cell=[u1,u2,u3],pbc=boundary)
        
        f.close()
        return atoms

    def imd_2_ase(fileName='IMD_conf_out', boundary = [1,1,1]):
        
        f = open(fileName,'r')
        line = f.readline()
        data = elements_in_line(line)
        while (data[0] != '#E'):
            print data[0]
            if (data[0] == '#X'):
                u1 = float(data[1])
            if (data[0] == '#Y'):
                u2 = float(data[2])
            if (data[0] == '#Z'):
                u3 = float(data[3])
            line = f.readline()
            data = elements_in_line(line)
        
        indices = []
        types = []
        mass = []
        pos = []
        
        try:
            while 1:
                line = f.readline()
                data = elements_in_line(line)
                indices.append(int(data[0])-1)
                types.append(int(data[1]))
                mass.append(float(data[2]))
                pos.append([float(data[3]), float(data[4]), float(data[5])])
        except IndexError:
            pass
        
        atomlist = []
        
        for i in range(len(indices)):
            i_index = indices.index(i)
            #        print pos[i_index]
            #        raw_input()
            atomlist.append(ase.Atom(position=pos[i_index],symbol=int(types[i_index])))
        
        atoms = ase.Atoms(atomlist,cell=[u1,u2,u3],pbc=boundary)
        
        f.close()
        return atoms

            
    def ase_2_con(atoms, filter=None, fileName=None, sort=False, pictures=1):
        if sort:
            sortTypes(atoms)
        
        pos = atoms.get_positions()
        cell = atoms.get_cell()
        types = atoms.get_atomic_numbers()
        
        ## all atoms are free
        if filter == None:
            filter = zeros(len(types))
        else:
            filter = abs(filter-1) #have to invert filter
        
        if fileName == None:
            fileName = raw_input("Enter file name (.con is added): ")
        fileName = fileName+".con"
        f = open(fileName,'w')
        
        locked = sum(filter)
        free = len(filter)-locked
        
        types_con = []
        types_nr_con = []
        last_type = None
        
        for type in types:
            if type != last_type:
                types_con.append(type)
                types_nr_con.append(1)
                last_type = type
            else:
                types_nr_con[-1] = types_nr_con[-1]+1
        
        f.write("10000 RANDOM NUMBER SEED\n") # for good random numbers
        f.write("0.0000 TIME\n") # starting time for simulation
        f.write(("%18f  %18f   %18f\n") % (cell[0][0], cell[1][1], cell[2][2])) # box size
        f.write(("%18f  %18f   %18f\n") % (90., 90., 90.)) # angles unitcell
        f.write("0 0\n") # pressure stuff
        f.write(str(locked)+ "  "+str(free)+" "+str(pictures)+"\n")
        f.write("  "+str(int(len(types_con)))+"\n") # number of different kinds of spicies
        
        for nr in types_nr_con:
            f.write(str(int(nr))+" ") # number of the different spicies
        f.write("\n")
        
        for temp in range(0,len(types_con)):
            f.write("1.0 ") # mass diffrent atom types
        f.write("\n")
        
        i = 0
        for type in range(0, len(types_con)):
            f.write(atom.elements[types_con[type]]['symbol']+"\n") # element types
            f.write("Coordinates of Component "+str(type+1)+"\n")
            for temp in range(0,types_nr_con[type]):
                f.write(("%18.15f  %18.15f   %18.15f   %i   %i\n") % (pos[i][0], pos[i][1], pos[i][2], filter[i], i+1))
                i = i+1
        return
            
            
    # convert atoms object
    def ase_2_eon(ase_atoms, free = None):
        
        e_atoms = et.Atoms(0)
        e_atoms.r = ase_atoms.get_positions()
        e_atoms.num_atoms = ase_atoms.get_number_of_atoms()
        e_atoms.box = ase_atoms.get_cell()
        e_atoms.free = free or np.ones(e_atoms.num_atoms,dtype='int')
        e_atoms.atomic_numbers = ase_atoms.get_atomic_numbers()
        e_atoms.mass = ase_atoms.get_masses()
        
        return e_atoms

    def ase_anim_xyz(confs, filename='Anim.xyz'):
        f_a = open(filename,'w')
        
        ## needed to make an xyz animation
        symbol = '/'
        stars = ''
        new_star = -1

        iteration = 1
        for conf in confs:
            atom_pos = []
            atomsType = conf.get_atomic_numbers()
            posAll = conf.get_positions()
            for posOne in posAll:
                atom_pos.append("%8.5f   %8.5f   %8.5f"%(posOne[0], posOne[1], posOne[2]))

            ## animation file
            if symbol == '/':
                symbol = '-'
                if new_star == 1:
                    stars = stars+'*'
                    new_star = 0
                else:
                    new_star += 1
                
            elif symbol == '-':
                symbol = '\\'
            elif symbol == '\\':
                symbol = '|'
            elif symbol == '|':
                symbol = '/'
            f_a.write(str(len(conf))+'\n')
            f_a.write(symbol+stars+'\n')

            i = 0
            for posOne in atom_pos:        
                f_a.write('   '+atom.elements[i]['symbol']+'   '+posOne+'\n')
                i = i+1

        f_a.close()
        return

    def ase_2_anim_xyz_interpolate(atoms1, atoms2, frames=15):
        initPos = atoms1.get_positions()
        movementStep = diff(atoms1, atoms2)[0]/float(frames)
        confs = []
        for i in range(0,frames):
            atomTemp = Atoms(atoms1)
            atomTemp.set_positions(initPos+i*movementStep)
            correctPos(atomTemp)
            confs.append(atomTemp)

        ase_2_anim_xyz(confs)
        return

    def ase_2_lammps(atoms, fileName=None):    
        pos = atoms.get_positions()
        cell = atoms.get_cell()
        types = atoms.get_atomic_numbers()

        print "write",pos[0]

        if fileName == None:
            fileName = raw_input("Enter file name (.conf is added): ")
        fileName = fileName+".conf"
        f = open(fileName,'w')

        types_con = []
        types_nr_con = []
        last_type = None
        types_masses = []
        for type in types:
            print type
            if type != last_type:
                types_con.append(type)
                types_nr_con.append(1)
                last_type = type
                types_masses.append(1.0)
            else:
                types_nr_con[-1] = types_nr_con[-1]+1

        f.write("Name\n\n")
        f.write(("%i  atoms\n") % len(pos))

        f.write(("%i  atom types\n") % len(types_masses))

        f.write(("0.0 %18f  xlo xhi\n") % cell[0][0]) # box size
        f.write(("0.0 %18f  ylo yhi\n") % cell[1][1]) # box size
        f.write(("0.0 %18f  zlo zhi\n") % cell[2][2]) # box size
        
        
        f.write(("\nMasses\n\n"))
        for nr in range(len(types_masses)):
            f.write(("%i %i\n") % (nr+1, types_masses[nr])) # number of the different spicies
        
        f.write(("\nAtoms\n\n"))
        i = 0
        for atomType in range(0, len(types_con)):
            for temp in range(0,types_nr_con[atomType]):
                f.write(("%4i %i %18.15f  %18.15f   %18.15f\n") % (i+1, atomType+1,
                                                                   pos[i][0], pos[i][1], pos[i][2]))
                i = i + 1
        return

    def ase_2_imd(atoms, filename='IMD_conf_out'):
        cell = atoms.get_cell()
        atomsType = atoms.get_atomic_numbers()
        posAll = atoms.get_positions()
        
        f = open(filename,'w')
        
        f.write(("#X %18f %18f %18f\n")%(cell[0][0], 0, 0)) # box size
        f.write(("#Y %18f %18f %18f\n")%(0, cell[1][1], 0)) # box size
        f.write(("#Z %18f %18f %18f\n")%(0, 0, cell[2][2])) # box size        
        
        f.write("#E\n")  
        
        i = 0
        for posOne in posAll:
            f.write("%4i %2i %8.5f %8.5f   %8.5f   %8.5f\n" % (i+1, atomsType[i], 4.0, posOne[0], posOne[1], posOne[2]))
            i = i + 1
        f.close()
        return

    def ice_xyz_2_ase(fileName, boundary = [1,1,1]):
        
        f = open(fileName,'r')
        line = f.readline()
        nrMolecules = int(elements_in_line(line)[0])
        
        line = f.readline()
        boxData = elements_in_line(line)
        
        u1 = float(boxData[0])
        u2 = float(boxData[1])
        u3 = float(boxData[2])
        
        atoms_read = []
        
        for i in range(nrMolecules):
            line = f.readline()
            data = elements_in_line(line)
            atoms_read.append([data[0], [float(data[1]), float(data[2]), float(data[3])]])
        
        atomlist = []
        
        for i in range(len(atoms_read)):
            atomlist.append(ase.Atom( symbol=atoms_read[i][0], position=atoms_read[i][1] ))
        
        atoms = ase.Atoms(atomlist,cell=[u1,u2,u3],pbc=boundary)
        f.close()
        return atoms

## ASE TO AND FROM SOMETHING END

##-----------------------
## CON TO SOMETHING START

#def con_2_ase(fileName = None, boundary = [1,1,1]):
#    if fileName == None:
#        fileName = raw_input("Enter file name (.con is added): ")
#        fileName = fileName+".con"
#
#    f = open(fileName,'r')
#	
#    ## seed to random
#    line = f.readline()
#
#    ## time stuff
#    line = f.readline()
#
#    line = f.readline()
#    basis_vec = elements_in_line(line)
#
#    line = f.readline()
#    angle_basis_vec = elements_in_line(line)
#    if int(float(angle_basis_vec[0])) != 90 or int(float(angle_basis_vec[1])) != 90 or int(float(angle_basis_vec[2])) != 90:
#        print "Only orthogonal lattice vectores are accepted!"
#        atoms = None
#	
#    else:
#        ## for storing the atoms
#        atomlist = Atoms([])
#        filter = []
#        ## pressure stuff
#        line = f.readline()
#        ## stuff
#        line = f.readline()
#		## number of different elements
#        line = f.readline()
#        nrTypes = int(elements_in_line(line)[0])
#        ## number of atoms of the different kinds
#        line = f.readline()
#        nrAtomsPerType = elements_in_line(line)
#        ## the different types mass
#        line = f.readline()
#        for i in range(0,nrTypes):
#            ## type of atom
#            line = f.readline()
#            type = elements_in_line(line)[0]
#            atomType = atom.elements[type]['number']
#            
#            ## text telling the components number
#            line = f.readline()
#            for j in range(0,int(nrAtomsPerType[i])):		
#                ## info about a certain atom
#                line = f.readline()                
#                atomInfo = elements_in_line(line)
#                atomlist.append(Atom(symbol=atomType,position=(float(atomInfo[0]),float(atomInfo[1]),float(atomInfo[2]))))
#                filter.append(int(atomInfo[3]))
#    atoms=Atoms(atomlist,cell=[float(basis_vec[0]),float(basis_vec[1]),float(basis_vec[2])],pbc=boundary)
#    f.close()
#    return (atoms, filter)
##return (atoms, abs(array(filter)-1))



## CON TO SOMETHING END
##-----------------------
## SOMETHING TO ASAP START

#def contcar2Asap(fileName = 'CONTCAR', boundary = [1,1,1]):
#    
#    fileName = fileName
#    f = open(fileName,'r')
#
#    ## atom type
#    line = f.readline()
#    types = elements_in_line(line)
#    
#    atomTypes = zeros(len(types))
#    index = 0
#    for type in types:
#        atomTypes[index] = atom.elements[type]['number']
#        index = index+1
#        
#    ## box size
#    line = f.readline()
#    boxSize = float(elements_in_line(line)[0])
#
#    ## unit vector 1
#    line = f.readline()
#    data = elements_in_line(line)
#    u1 = array([float(data[0]), float(data[1]), float(data[2])])*boxSize
#    ## unit vector 2
#    line = f.readline()
#    data = elements_in_line(line)
#    u2 = array([float(data[0]), float(data[1]), float(data[2])])*boxSize
#    ## unit vector 3
#    line = f.readline()
#    data = elements_in_line(line)
#    u3 = array([float(data[0]), float(data[1]), float(data[2])])*boxSize
#
#    ## number of atoms
#    line = f.readline()
#    nrs = elements_in_line(line)
#
#    atomNrs = zeros(len(types))
#    index = 0
#    for nr in nrs:
#        atomNrs[index] = int(nr)
#        index = index+1
#
#    ## crap that does not matter
#    line = f.readline()
#    ## if scaled or cartesian positions are used
#    line = f.readline()
#    coorType = elements_in_line(line)
#    if coorType[0] == 'Cartesian':
#        s1 = array([1., 0, 0])
#        s2 = array([0, 1., 0])
#        s3 = array([0, 0, 1.])
#    else:
#        s1 = u1
#        s2 = u2
#        s3 = u3
#        
#    ## for storing the atoms
#    atomlist = Atoms([])
#    filter = []
#
#    for index in range(0,len(atomNrs)):
#        for i in range(0,atomNrs[index]):
#        ## info about a certain atom
#            line = f.readline()                
#            atomInfo = elements_in_line(line)
#            pos = array(float(atomInfo[0])*s1+float(atomInfo[1])*s2+float(atomInfo[2])*s3)
#            atomlist.append(Atom(position=pos,symbol=int(atomTypes[index])))
#            if atomInfo[3] == 'T' and atomInfo[4] == 'T' and atomInfo[5] == 'T':
#                filter.append(1)
#            else:
#                filter.append(0)
#
#    atoms = Atoms(atomlist,cell=[u1,u2,u3],pbc=boundary)
#    
#    f.close()
#    return (atoms, filter)
#
#def poscar2Asap(fileName = 'POSCAR', boundary = [1,1,1]):
#    return contcar2Asap(fileName,boundary)
#
#def contcarEnergy2Asap(fileName = 'POSCAR', boundary = [1,1,1]):
#    energy = 0
#
#    (atoms, filter) = contcar2Asap(fileName,boundary)
#
#    fileName = fileName
#    f = open(fileName,'r')
#    
#    lines = f.readlines()
#    for line in lines:
#        data = elements_in_line(line)
#        if len(data) > 0:
#            if data[1:3] == ['Relative','Energy']:
#                energy = float(data[0])
#                break
#            
#    f.close()
#    return (atoms, filter, energy)



#def xyz2Asap(fileName, boundary = [1,1,1]):
#    f = open(fileName,'r')
#    line = f.readline()
#    data = elements_in_line(line)
#    atomNrs = int(data[0])
#    
#    atoms_read = {}
#    
#    line = f.readline()
#    for i in range(atomNrs):
#        line = f.readline()
#        data = elements_in_line(line)
#        
#        try:
#            atomNr = int(data[0])
#        except ValueError:
#            atomNr = atom.elements[data[0]]['number']
#
#        atoms_read[i] = [atomNr, [float(data[1]), float(data[2]), float(data[3])]]
#    
#    atomlist = []
#    
#    keys = atoms_read.keys()
#    keys.sort()
#    for key in keys:
#        atomlist.append(Atom( symbol=atoms_read[key][0], position=atoms_read[key][1] ))
#    
#    atoms = Atoms(atomlist,cell=[1000,1000,1000],pbc=boundary)
#    
#    f.close()
#    return atoms
