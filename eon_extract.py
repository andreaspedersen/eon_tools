#!/usr/bin/python
import eon_read
import atom
import numpy
import molecule
import math
from operator import itemgetter

######################################################################################################

def extract_traj_atom(index):
    states = eon_read.read_statelist()
    dynamics = eon_read.read_dynamics()

    traj_atom = []
    for d in dynamics:
        pos = states[d['reactantid']]['reactant_config'].r
        traj_atom.append(pos[index])
    return traj_atom

######################################################################################################

def extract_atom_position_all_states(index,statelist=None):
    if statelist ==None:
        statelist = eon_read.read_statelist(read_reactant=True)
    
    positions = len(statelist)*[None]
    for state in statelist:
        positions [ state['stateid'] ] = extract_atom_position(index, state['reactant_config'] )
    
    return positions

######################################################################################################

def extract_atom_position(index,configuration):
    
    if not isinstance(index,list):
        atomlist = [ index ]
    else:
        atomlist = index
    
    box = configuration.box
    r = numpy.zeros((3))
    
    for atomid in atomlist:
        r += _fit_box( configuration.r[ atomid ] - configuration.r[ atomlist[0] ] , box )
    
    r = r / len(atomlist)
    r += configuration.r[ atomlist[0] ]

    for i in range(3):
        while r[i] > box[i][i]: r[i] -= box[i][i]
        while r[i] < 0.0: r[i]+=box[i][i]

    return r

######################################################################################################

def extract_atom_displacement(index,configuration1,configuration2):
    '''
    Returns the displacement of the center of the id's listed in the index list
    going from configuration1 to configuration2
    '''
    
    pos1 = extract_atom_position(index,configuration1)
    pos2 = extract_atom_position(index,configuration2)
    
    displacement = _fit_box( (  pos2 - pos1 ) , configuration1.box )
    
    return displacement

######################################################################################################

def extract_com_displacement(index,configuration1,configuration2):
    '''
    Returns the displacement of the COM of the atoms listed in index going from
    configuration1 to configuration2
    '''
    
    com1 = extract_com(index,configuration1)
    com2 = extract_com(index,configuration2)
    
    displacement = _fit_box( ( com2 - com1 ), configuration1.box)
    
    return displacement

######################################################################################################

def extract_displacement_atom(index):
    state = eon_read.read_state(0,read_configurations=True)
    box = state['reactant_config'].box
    traj_atom = extract_traj_atom(index)
    
    diff_all = []
    for i in range(len(traj_atom)-1):
        diff = traj_atom[i+1]-traj_atom[i]
        diff_all.append(_fit_box(diff, box))
    return diff_all

######################################################################################################

def extract_energies_path():
    dynamics = eon_read.read_dynamics()
    states = eon_read.read_statelist()
    try :
        sb = eon_read.read_superbasins()
    except AssertionError:
        sb = None

    energies = []
    
    for i in range(len(dynamics)):
        # Min energy
        if dynamics[i]['steptype'] == 'kmc':
            energies.append(states[dynamics[i]['reactantid']]['energy'])
        
        elif dynamics[i]['steptype'] == 'sb':
            energies.append(sb[dynamics[i]['sbid']]['min_energy'])
        
        else:
            print "Unknown dynamics step at "+str(i)
            stop

    return energies
"""
def extract_energies_path():
    dynamics = eon_read.read_dynamics()
    states = eon_read.read_statelist()
    sb = eon_read.read_superbasins()
    energies = []

    for i in range(len(dynamics)):
        try:
            energies.append(sb[dynamics[i]['superbasin_id']]['min_energy'])
        except KeyError:
            energies.append(states[dynamics[i]['reactantid']]['energy'])
    return energies
"""
######################################################################################################

def extract_com(index,configuration):
    '''
    Returns the COM displacement of the list of atoms in index.
    '''
    try:
        natoms = len(index)
    except:
        natoms = 1
        index = [ int(index) ]
    
    origin = configuration.r[ index[0] ]
    com = numpy.zeros(3)
    #calculate everything w.r.t atom index[0]
    totalmass = sum( [ configuration.mass[ i ] for i in index] )
    for i in range( len(index) -1 ):
        atid=index[i+1]
        com += configuration.mass[atid] * _fit_box( ( configuration.r[atid] - configuration.r[ index[0] ] ), configuration.box )
    com = com*1.0/totalmass
    com +=origin

    for i in range(3):
        if (com[i] > configuration.box[i][i]):
            com[i] -= configuration.box[i][i]
        elif (com[i] < 0.0):
            com[i] += configuration.box[i][i]
    return com

######################################################################################################

def extract_energies_path_with_saddles():
    dynamics = eon_read.read_dynamics()
    states = eon_read.read_statelist()
    
    try :
        sb = eon_read.read_superbasins()
    except AssertionError:
        sb = None

    energies = []
    for i in range(len(dynamics)):
        # Min energy
        if dynamics[i]['steptype'] == 'kmc':
            energies.append(states[dynamics[i]['reactantid']]['energy'])
        
        elif dynamics[i]['steptype'] == 'sb':
            energies.append(sb[dynamics[i]['sbid']]['min_energy'])
        
        else:
            print "Unknown dynamics step at "+str(i)
            stop
    
        # SP energy
        energies.append(eon_read.read_process(dynamics[i]['reactantid'],dynamics[i]['processid'])['saddle_energy'])
    return energies
"""
def extract_energies_path_with_saddles():
    dynamics = eon_read.read_dynamics()
    states = eon_read.read_statelist()
    sb = eon_read.read_superbasins()

    energies = []
    for i in range(len(dynamics)):
        try:
            energies.append(sb[dynamics[i]['superbasin_id']]['min_energy'])
        except KeyError:
            energies.append(states[dynamics[i]['reactantid']]['energy'])

        energies.append(eon_read.read_process(dynamics[i]['reactantid'],dynamics[i]['processid'])['saddle_energy'])
    return energies
"""
######################################################################################################

def extract_barriers():
    states = eon_read.read_statelist()
    
    
    barriers = []
    for step in eon_read.dynamics():
        min_e = states[step['reactantid']]['energy']
        sp_e = eon_read.read_process(step['reactantid'],step['processid'])['saddle_energy']
        barriers.append(sp_e-min_e)
    return barriers

######################################################################################################
#Function extracts the first N nearest neighbours of molecule molid

def extract_molecule_nearest_neighbours(config,molid,n,radius = 6.0):
    eon_read.read_molecules(config,read_com=True)
    config.create_verletList(radius)
    for mol in config.molecules:
        config.verletList.append(mol._molid,mol._com)
    
    neighbourdict = { }
    mol = config.molecules[molid]
    for neighbourid in config.verletList.neighbours(molid):
        nmol = config.molecules[neighbourid]
        if nmol._molid == molid:
            continue
        distance = numpy.linalg.norm(_fit_box(nmol._com - mol._com, config.box) )
        neighbourdict[nmol._molid] = distance

    #sort by value
    returnlist = sorted(neighbourdict.items(), key=itemgetter(1))
    if len(returnlist) < n and radius < max( [ config.box[i][i]/3 for i in range(3) ]):
        returnlist = extract_nearest_neighbours(config,molid,n,radius*2)

    if len(returnlist) < n:
        sys.stderr.write("WARNING: only %d neighbours could be found for mol %d (%d requested)\n" % len(returnlist), molid, n)
    return returnlist[:n]


######################################################################################################

def moved_atoms(reactant, product, cutoff):
    assert( len(product)==len(reactant) )
    moved_atom_ids = []
    
    for iatom in range(len(reactant)):
        r = numpy.linalg.norm( _fit_box( reactant.r[iatom] - product.r[iatom], reactant.box) )
        if r > cutoff:
            moved_atom_ids.append(iatom)

    return moved_atom_ids

######################################################################################################

def moved_molecules(reactant, product, cutoff):
    try:
        reactant.molecules
    except:
        eon_read.read_molecules(reactant)
    try:
        product.molecules
    except:
        eon_read.read_molecules(product)

    assert product.molecules._num_molecules == reactant.molecules._num_molecules

    moved_mol_ids = []
    for imol in range( reactant.molecules._num_molecules ):
        assert reactant.molecules._list[imol]._name == product.molecules._list[imol]._name
        for iatom in reactant.molecules._list[imol]._atomids:
            r = numpy.linalg.norm( reactant.r[iatom] - product.r[iatom] , reactant.box)
            if r > cutoff:
                moved_mol_ids.append(imol)
                break


######################################################################################################

def molecule_movement(reactant, product, disp_cutoff=0.0, angle_cutoff=0.0):
    
    #determine molecules
    eon_read.read_molecules(reactant,read_com=True,read_orientation=True)
    eon_read.read_molecules(product,read_com=True,read_orientation=True)
    
    returnlist = []
    for molid in range(reactant.molecules._num_molecules):
        if not reactant.molecules._list[molid]._free:
            continue
        movementdict = {'molid': molid}
        #calculate COM displacement
        disp = numpy.linalg.norm(_fit_box(reactant.molecules._list[molid]._com - product.molecules._list[molid]._com, reactant.box ) )
        movementdict['disp'] = disp
    
        #calculate angle of rotation
        n1 = reactant.molecules._list[molid]._orientation
        n2 =  product.molecules._list[molid]._orientation
        dotproduct = n1.dot(n2)
        crossproduct = numpy.linalg.norm( numpy.cross(n1,n2) )
        angle = math.atan2( crossproduct,dotproduct)*180.0/math.pi
        movementdict['angle'] = angle
    
        #if water molecule there is one more angle
        if reactant.molecules[molid]._name=="H2O":
            n1 = reactant.molecules[molid]._orientation2
            n2 = product.molecules[molid]._orientation2
            dotproduct = n1.dot(n2)
            crossproduct = numpy.linalg.norm( numpy.cross(n1,n2) )
            angle = math.atan2( crossproduct,dotproduct)*180.0/math.pi
            if angle > movementdict['angle']:
                movementdict['angle'] = angle
    
        #add to return list
        if disp >= disp_cutoff or math.fabs(angle) >= angle_cutoff:
            returnlist.append(movementdict)

    #return
    return returnlist

######################################################################################################

def extract_effective_barriers():
    dynamics = eon_read.read_dynamics()
    states = eon_read.read_statelist()

    try :
        sb = eon_read.read_superbasins()
    except AssertionError:
        sb = None
    
    barriers = []
    for i in range(len(dynamics)):
        # Min energy
        if dynamics[i]['steptype'] == 'kmc':
            min_e = states[dynamics[i]['reactantid']]['energy']
        
        elif dynamics[i]['steptype'] == 'sb':
            min_e = sb[dynamics[i]['sbid']]['min_energy']
        
        else:
            print "Unknown dynamics step at "+str(i)
            stop
        
        # SP energy
        sp_e = eon_read.read_process(dynamics[i]['reactantid'],dynamics[i]['processid'])['saddle_energy']

        # Barrier
        barriers.append(sp_e-min_e)
    return barriers
"""
def extract_effective_barriers():

    states = eon_read.read_statelist()
    sb = eon_read.read_superbasins()
    print "read_superbasins"
    barriers = []
    #for i in range(len(dynamics)):
    for step in eon_read.dynamics():
        if step['steptype']=='sb':
            min_e = sb[ step['sbid']]['min_energy']
        else:
            min_e = states[step['reactantid']]['energy']
        sp_e = eon_read.read_process(step['reactantid'],step['processid'])['saddle_energy']
        barriers.append(sp_e-min_e)
    return barriers
"""
######################################################################################################

def extract_time():
    dynamics = eon_read.read_dynamics()

    time = []
    for i in range(len(dynamics)):
        time.append(dynamics[i]['totaltime'])
    return time

######################################################################################################

def extract_time_spend_in_states():
    dynamics = eon_read.read_dynamics()
    
    try :
        sb = eon_read.read_superbasins()
    except AssertionError:
        sb = None

    end_state = 0
    t_last_step = 0
    
    state_time = {}
    
    for i in range(len(dynamics)):
        reactant_id = dynamics[i]['reactantid']
        
        t_step = dynamics[i]['totaltime'] - t_last_step
        t_last_step = dynamics[i]['totaltime']

        if reactant_id > end_state:
            end_state = reactant_id

        if dynamics[i]['steptype'] == 'sb':
            reactant_id = sb[dynamics[i]['sbid']]['min_state_id']

        try:
            state_time[reactant_id] = state_time[reactant_id] + t_step
        except:
            state_time[reactant_id] =  t_step

    for i in range(end_state):
        try:
            state_time[i]
        except:
            state_time[i] = 0.

    return state_time

######################################################################################################

def extract_visit_to_states():
    dynamics = eon_read.read_dynamics()
    
    try :
        sb = eon_read.read_superbasins()
    except AssertionError:
        sb = None
    
    end_state = 0
    
    state_visits = {}
    
    for i in range(len(dynamics)):
        reactant_id = dynamics[i]['reactantid']

        if reactant_id > end_state:
            end_state = reactant_id
        
        if dynamics[i]['steptype'] == 'sb':
            reactant_id = sb[dynamics[i]['sbid']]['min_state_id']
        
        try:
            state_visits[reactant_id] = state_visits[reactant_id] + 1
        except:
            state_visits[reactant_id] =  1
    
    for i in range(end_state):
        try:
            state_visits[i]
        except:
            state_visits[i] = 0.
    
    return state_visits

######################################################################################################

def extract_state_energies(states=None):
    e = []
    if states == None:
        f = open('state_table')
        lines = f.readlines()
        for line in lines:
            e.append(float(line.split()[1]))
    else:
        for state in states:
            e.append(state['energy'])

    return e

######################################################################################################

def _fit_box(diff, box):
    for i in range(3):
        if abs(diff[i]) > box[i][i]/2.:
            if diff[i] > box[i][i]/2.:
                diff[i] = diff[i]-box[i][i]
            else:
                diff[i] = diff[i]+box[i][i]

            if box[i][i] * 0.45 < diff[i]:
                print 'Warning displace close to half a box length!'
                
    return diff
