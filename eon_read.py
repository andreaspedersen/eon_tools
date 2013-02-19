#!/usr/bin/python
import os
import sys
import pickle

import os, sys
import numpy as np
import ConfigParser

import atom
import molecule
import eon_extract

class AKMC_proc_data(object):
    """
    Takes a path to data created by eOn AKMC calculation and makes the data simply
    accessable through functioncalls.
    """

    def __init__(self, path):
        self.path = path
        self.nrStates = self._findNumStates()
        self.nrProcesses = self._findNumProc()
        self.totalNrProcesses = self._findTotProc()
        
        
    def _findNumStates(self):
        fName = 'state_table'
        path = os.path.join(self.path, 'states')
        try:
            res = len(np.loadtxt(os.path.join(path, fName), usecols=(0,)))
        except IOError:
            print 'file state_table is empty. Assuming single state.'
            res = 1
        return res
        
    def _findNumProc(self):
        procs = np.empty(self.nrStates, dtype='int')
        fName = 'processtable'
        for state in xrange(self.nrStates):
            path = os.path.join(self.path, 'states', str(state))
            out = np.loadtxt(os.path.join(path, fName),
                             skiprows=1, usecols=(0,))
            procs[state] = out.size
        return procs
        
    def _findTotProc(self):
        result = 0
        for key in xrange(self.nrStates):
            result += self.nrProcesses[key]
        return result
     
    def _getProcessTable(self, state):
        if state < self.nrStates:
            path = os.path.join(self.path, 'states', str(state), 'processtable')
            return path
            
    
    def _getOrderBy(self, f, index):
        """
        Takes column in a file (proctable) and returns sorted indexed
        """
        arr = np.loadtxt(f,skiprows=1,usecols=(index,))
        return np.argsort(arr)
        
    def getProcIndexByLowestSaddleEnergy(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = self._getOrderBy(f, 1)
        f.close()
        return res
            
    def getProcIndexByLowestProductEnergy(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = self._getOrderBy(f, 4)
        f.close()
        return res
            
    def getSaddleEnergys(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(1,))
        f.close()
        return res
            
    def getProductEnergys(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(4,))
        f.close()
        return res
        
    def getProductPicked(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(3,))
        f.close()
        return res
        
    def getBarrier(self, state):
        f = open(self._getProcessTable(state), 'r')
        res = np.loadtxt(f,skiprows=1,usecols=(6,))
        f.close()
        return res
            
    def getReactantEnergy(self, state):
        filePath = os.path.join(self.path, 'states', str(state), 'info')
        f = open(filePath, 'r')
        cfg = ConfigParser.ConfigParser()
        cfg.readfp(f)
        f.close()
        return cfg.getfloat('MetaData', 'reactant energy')
            
            
    def getLowestProductEnergyOfAll(self, nrProc):
        energy = np.zeros(0)
        states = np.zeros(0, dtype='int')
        process = np.zeros(0, dtype='int')
        for state in xrange(self.nrStates):
            f = open(self._getProcessTable(state), 'r')
            arr = np.loadtxt(f,skiprows=1,usecols=(4,))
            f.close()
            energy = np.hstack((energy,arr))
            arg = np.argsort(energy)
            states = np.hstack((states,np.ones(arr.size,dtype='int')*state))
            process = np.hstack((process,
                                 np.array(range(arr.size),dtype='int')))
            energy = energy[arg][:nrProc]
            states = states[arg][:nrProc]
            process = process[arg][:nrProc]
        return energy, states, process
        
    def getReac(self, state, proc):
        if state < self.nrStates and proc < self.nrProcesses[state]:
            fName = 'reactant.con'
            path = os.path.join(self.path, 'states', str(state))
            return os.path.join(path, fName)
        print "State or process doesn't exist"
    
    def getProd(self, state, proc):
        """
        Get all products for specific state
        """
        if state < self.nrStates and proc < self.nrProcesses[state]:
            fName = 'product_%i.con' %(proc,)
            path = os.path.join(self.path, 'states', str(state), 'procdata')
            return os.path.join(path, fName)
        print "State or process doesn't exist"
    
    def getReacProd(self, state, proc):
        """
        Get all processes for specific state
        """

        return self.getReac(state,proc), self.getProd(state,proc)
    
    def getAllReacProd(self):
        """
        Returns all processes in all proctables
        """
        for state in xrange(self.nrStates):
            for proc in xrange(self.nrProcesses[state]):
                yield self.getReacProd(state, proc), (state, proc)
                
    def getInitialReac(self, state):
        if state < self.nrStates:
            fName = 'reactant.con'
            path = os.path.join(self.path, 'states', str(state))
            return os.path.join(path, fName)
        print "State doesn't exist"
                
    def stateProcToIndex(self, state, proc):
        return self.nrProcesses[:state].sum() + proc
        
    def indexToStateProc(self, index):
        states = 0
        if index < 0:
            raise IndexError("Negative indexes not allowed.")
        for state in xrange(self.nrStates):
            
            if index < states + self.nrProcesses[state]:
                proc = index - states
                return state, proc
            states += self.nrProcesses[state]
        raise IndexError("Index %i is out of range. Highest value is %i" %(index, states))


######################################################################################################

def read_state_table(path='states'):
    '''
        Returns a dictionary containing state energies
        '''
    state_table = {}
    state_table_file = _open_file_reading(os.path.join(path,"state_table"))
    if state_table_file != None:
        lines = state_table_file.readlines()
        state_table_file.close()
        
        for line in lines:
            data = line.split()
            state_table[int(data[0])] = float(data[1])
    return state_table    

######################################################################################################

def read_proc_configs(stateid,procid,path='states'):
    '''
    Returns a dictionary containing atom objects of the reactant,saddle and procuct configurations of the specified process
    '''
    processpath = os.path.join(path,str(stateid),"procdata")
        
    reactantcon = _read_con(os.path.join(processpath,"reactant_"+str(procid)+".con"))
    saddlecon   = _read_con(os.path.join(processpath,"saddle_"+str(procid)+".con"))  
    productcon  = _read_con(os.path.join(processpath,"product_"+str(procid)+".con"))   
       
    return {'reactant_config': reactantcon,'saddle_config': saddlecon, 'product_config': productcon }


######################################################################################################

def read_process(stateid, procid, path='states', read_configurations=False, read_all=False):
    '''
    Function returns a dictionary containing all information about the given process.
    If read_configurations = True, atom objects will be included containing the configurations of the reactant,saddle and product
    '''
    statepath = os.path.join( path,str(stateid) )
    if not os.path.isdir(path):
        print "error: %s is not a valid statepath" % path
        return
    
    
    processtable = _open_file_reading(os.path.join(path,str(stateid),"processtable"))
    barrier = None
    procdict = { }
    if processtable !=None:
        lines = processtable.readlines()
        processtable.close()
        assert(procid<len(lines))
        saddle_energy = float(lines[procid+1].split()[1] )
        prefactor = float(lines[procid+1].split()[2] )
        productid = int(lines[procid+1].split()[3] )        
        product_energy = float(lines[procid+1].split()[4] )
        product_prefactor = float(lines[procid+1].split()[5] )
        barrier =float( lines[procid+1].split()[6]  )
        rate = float( lines[procid+1].split()[7]  )
        repeats = int( lines[procid+1].split()[8]  )


        procdict = {'procid' : int(procid),
                    'rate'   : rate,
                    'repeats' : repeats, 
                    'barrier': barrier,
                    'prefactor': prefactor,
                    'productid': productid,
                    'reactantid': stateid,
                    'saddle_energy': saddle_energy,
                    'product_energy': product_energy,
                    'reactant_energy': ( saddle_energy - barrier ),
                    'product_prefactor': product_prefactor
                    }
        if read_configurations == True:
            configurations_dict = read_proc_configs(stateid,procid,path)
            procdict = dict(procdict.items() + configurations_dict.items() )

        if read_all == True:
            procdict['modefile']    = _read_proc_modefile(path,stateid,procid)
            procdict['resultsfile'] = _read_proc_resultsfile(path,stateid,procid)
    return procdict
    

######################################################################################################
    
def read_state(stateid,path='states',read_reactant=True, read_processes=False, read_configurations=False, read_all=False):
    '''
	Reads all data from a specified state into a dictionary and returns it.
    '''
    #print "reading state %d" % stateid
    sys.stderr.write("Reading state %d\n" % stateid)

    #read info file:
    statedict = _read_state_info(path,stateid)
        
    statedict['stateid'] = stateid
    statedict['numprocs'] = _determine_num_procs(path,stateid)

    if read_processes == True:
        statedict['proclist'] = [ ]
        for procid in range(statedict['numprocs'] ):
            statedict['proclist'].append ( read_process(stateid,procid,path,read_configurations,read_all) )
            i = 0
            for proc in statedict['proclist']:
                assert(proc['procid'] == i)
                i += 1

    if read_reactant == True:
        statedict['reactant_config'] = _read_con( os.path.join( path, str(stateid), "reactant.con" ) )
    
    if read_all == True:
        statedict['infofile'] = _read_infofile(path,stateid)
        statedict['search_results'] = _read_search_results(path,stateid)
    return statedict

######################################################################################################
def read_statelist(path='.', read_reactant=True, read_processes=False, read_configurations=False, read_all=False, read_configs_unknown=False, dump=0):
    '''
	Read data from all states and return it in a list containing dictionaries with information
    '''
    statelist = _attempt_read_pickle("states")
    
    if statelist is None:
        statespath = os.path.join(path,"states")
        
        if os.path.isdir(statespath)==False:
            print "ERROR: %s is not a valid eon path (is does not contain a states subdir" % path
            return

        if read_all==True:
            read_reactant=True
            read_processes=True
            read_configurations=True

        numstates = _determine_num_states(statespath)
        
        statelist = [ ]
        for stateid in range(numstates):
            statelist.append( read_state(stateid,statespath,read_reactant,read_processes,read_configurations,read_all) )

        if read_configs_unknown==True and read_configurations==False and read_reactant==True and read_processes==True:
        
            for state in statelist:
                sys.stderr.write("Reading unkown configs for state %d\n" % state['stateid'])
                for proc in state['proclist']:
                    proc['reactant_config'] = state['reactant_config']
                    if proc['productid']==-1:
                        processpath = os.path.join(statespath,str(state['stateid']),"procdata")
                        proc['product_config'] =  _read_con(os.path.join(processpath,"product_"+str(proc['procid'])+".con"))
                    else:
                        proc['product_config'] = statelist[ proc['productid'] ]['reactant_config']

        if dump:
            _dump_pickle("states", statelist)

    return statelist

######################################################################################################

def read_molecules(config,read_com= False, read_orientation = False):
    try:
        config.molecules
    except:
        config.molecules = molecule.MolList(config)
    
    if read_com==True:
        read_coms(config)
    if read_orientation == True:
        read_orientations(config)

    return

######################################################################################################

def read_coms(config):
    read_molecules(config)
    for mol in config.molecules:
        mol._com = eon_extract.extract_com(mol._atomids,config)

######################################################################################################

def read_orientations(config):
    read_molecules(config)
    for mol in config.molecules:
        mol.calc_orientation(config)

######################################################################################################

def read_dynamics(path='.',only_stateids=False, dump=0):
    
    dynamics = _attempt_read_pickle("dynamics")
    if dynamics is None:    
        dynamics = [ ]
        for step in dynamics(path,only_stateids):
            dynamics.append(step)

    if dump:
        _dump_pickle("dynamics", dynamics)
    return dynamics
"""
def read_dynamics(path='.',only_stateids=False, dump=0):

    dynamics = _attempt_read_pickle("dynamics")

    if dynamics is None:
        dynamics = [ ]
        #open dynamicsfile
        print path
        dynamicsfile = _open_file_reading(os.path.join(path,"dynamics.txt") )
        dynamicsfile.readline()
        dynamicsfile.readline()
        for line in dynamicsfile:
            splitline = line.split()

            #have not yet reached the end of file
            if len(splitline):
                stepid = int(splitline[0] )
                stepdict = {  }
                stepdict['reactantid'] = int(splitline[1])

                if only_stateids==False:

                    stepdict['processid'] = int(splitline[2])
                    stepdict['productid'] = int(splitline[3])
                    stepdict['timestep'] = float(splitline[4])
                    stepdict['totaltime'] = float(splitline[5])

                    #Superbasin step
                    if len(splitline) == 7:
                        stepdict['steptype'] = 'sb'
                        stepdict['sbid'] = int(splitline[6] )
                    elif len(splitline) == 8:
                        stepdict['steptype'] = 'kmc'
                    else:
                        print "Read error in dynamicsfile on line %s" % line

                dynamics.append(stepdict)

        if dump:
            _dump_pickle("dynamics", dynamics)

    return dynamics
"""
######################################################################################################

def dynamics(path='.',only_stateids=False):
    dynamicsfile = _open_file_reading(os.path.join(path,"dynamics.txt") )
    dynamics = [ ]
    dynamicsfile.readline()
    dynamicsfile.readline()
    for line in dynamicsfile:
        splitline = line.split()
        
        
        stepid = int(splitline[0] )
        stepdict = {  }
        stepdict['reactantid'] = int(splitline[1])
        
        if only_stateids==False:
            stepdict['stepid'] = stepid
            stepdict['processid'] = int(splitline[2])
            stepdict['productid'] = int(splitline[3])
            stepdict['timestep'] = float(splitline[4])
            stepdict['totaltime'] = float(splitline[5])
            
            #Superbasin step
            if len(splitline) == 7:
                stepdict['steptype'] = 'sb'
                stepdict['sbid'] = int(splitline[6] )
            elif len(splitline) == 8:
                stepdict['steptype'] = 'kmc'
            else:
                print "Read error in dynamicsfile on line %s" % line
    
        yield stepdict

    dynamicsfile.close()

######################################################################################################

def totaltime(path='.'):
    dynamicsfile = _open_file_reading(os.path.join(path,"dynamics.txt") )
    dynamicsfile.seek(0,2)
    filesize=dynamicsfile.tell()
    dynamicsfile.seek(max( filesize-2048 ,0 ) ,0 )
    lines=dynamicsfile.readlines()
    totaltime = float ( lines[-1].split()[5] )
    dynamicsfile.close()
    return totaltime

######################################################################################################

def read_superbasins(path='.', read_old_basins=True, dump=0):
    '''
	Reads superbasin information and returns a list of all superbasins.
	This list contains dictionaries with further information about the superbasins
    '''

    sblist = _attempt_read_pickle("superbasins")

    if sblist is None:
        superbasins = []
        superbasinpath = os.path.join(path,"superbasins")
        assert(os.path.isdir(superbasinpath))


        #determine the number of current superbasins
        for filename in os.listdir(superbasinpath):
            if filename == "storage":
                continue
            superbasinid=int(filename)
            sbdict = _read_superbasin(superbasinpath,filename)
            sbdict['id'] = superbasinid
            sbdict['is_active'] = True
            superbasins.append(sbdict)

        #determine all old superbasins
        oldsuperbasins = []
        if os.path.isdir( os.path.join(superbasinpath,"storage") ) and read_old_basins==True:
            for filename in os.listdir( os.path.join( superbasinpath,"storage") ):
                superbasinid = int(filename)
                sbdict = _read_superbasin( os.path.join( superbasinpath,"storage"), filename)
                sbdict['id'] = superbasinid
                sbdict['is_active'] = False
                oldsuperbasins.append(sbdict)
        
        #merge superbasins into list with sorted id's.
        maxid = max( [ i['id'] for i in superbasins+oldsuperbasins ] )
        sblist = [None] *( maxid + 1 )
        for sb in superbasins+oldsuperbasins:
            sblist[ int( sb['id'] ) ] = sb
            
        #determine minimum energy of each superbasin:
        statelist = read_statelist(path, read_reactant=False, read_processes=False, read_configurations=False)
        for sb in sblist:
            if sb == None:
                continue
            min_energy= 1e200
            for stateid in sb['states']:
                min_energy = min(min_energy, statelist[stateid]['energy'] )
                min_state_id = stateid
            sb['min_energy'] = min_energy
            sb['min_state_id'] = min_state_id

        if dump:
            _dump_pickle("superbasins", sblist)

            
    return sblist

######################################################################################################



















######################################################################################################
#   'PRIVATE' FUNCTIONS
######################################################################################################

def _attempt_read_pickle(filename):
    try:
        f = open(filename+".pickle",'r')
        data = pickle.load(f)
        f.close()
    except:
        data = None
    return data

def _dump_pickle(filename, data):
    f = open(filename+".pickle",'w')
    pickle.dump(data,f)
    return

def _read_superbasin(superbasinpath,sbfilename):
    sbfile = _open_file_reading(os.path.join(superbasinpath,sbfilename))

    line = sbfile.readline().split()
    numstates = len(line)
    states = [ ]
    for i in range(numstates):
        states.append( int(line[i]) )
        
    sbdict = { 'numstates': numstates, 'states':states}
    sbfile.close()    
    return sbdict     

######################################################################################################

def _determine_num_states(statepath):
    '''
	Determines the total number of states which the AKMC simulaton has visited. 
	This is simply determined from the number of lines in the state_table file.
    '''
    state_table = _open_file_reading(os.path.join(statepath,"state_table"))
    num_states = 0
    if state_table != None:
        lines = state_table.readlines()
        state_table.close()
        num_states = len(lines)
    return num_states

######################################################################################################

def _read_con(filename):
    config = atom.Atoms(0)
    config.read_con(filename)
    
    return config

######################################################################################################

def _determine_num_procs(statepath,stateid):
    '''
	Returns the number of unique processes leading out of a specified state.
    '''
    path = os.path.join(statepath,str(stateid) )
    if not os.path.isdir(path):
        print "error: %s is not a valid statepath" % path
    processtable = _open_file_reading(os.path.join(path,"processtable"))
    numprocs = 0
    if processtable != None:
        lines = processtable.readlines()
        processtable.close()
        numprocs = len(lines)-1
    return numprocs


######################################################################################################

def _read_state_info(statepath,stateid):    
    energy = 0
    good_saddles = -1
    bad_saddles = -1
    unique_saddles = -1
    
    try:
        infofile = _open_file_reading( os.path.join ( statepath, str(stateid), "info" ) )
    except IOError:
        print "Could not open file %s" %  os.path.join ( statepath, str(stateid), "info" )
        return
        
    for line in infofile:
        if len(line.split() ) > 2:
            if line.split()[0] == "reactant" and line.split()[1] == "energy":
                energy = float( line.split()[3] )
            elif line.split()[0] == "good_saddles":
                good_saddles = int( line.split()[2] )
            elif line.split()[0] == "bad_saddles":
                bad_saddles = int ( line.split()[2] )
            elif line.split()[0] == "unique_saddles":
                unique_saddles = int(line.split()[2])
    infofile.close()
    
    return {    'energy': energy,
                'unique_saddles' : unique_saddles,
                'bad_saddles': bad_saddles,
                'good_saddles': good_saddles
           }


######################################################################################################

def _read_infofile(statepath,stateid):
    infofile = _open_file_reading(os.path.join( statepath, str(stateid),"info") )
    lines = infofile.readlines()
    infofile.close()
    return lines


######################################################################################################

def _read_search_results(statepath,stateid):
    searchresultsfile = _open_file_reading(os.path.join(statepath, str(stateid), "search_results.txt"))
    lines=searchresultsfile.readlines()
    searchresultsfile.close()
    return lines


######################################################################################################

def _read_proc_modefile(statepath,stateid,procid):
    modefile = _open_file_reading( os.path.join(statepath,str(stateid),"procdata",str("mode_%d.dat" % procid) ) )
    lines = modefile.readlines()
    modefile.close()
    return lines

######################################################################################################

def _read_proc_resultsfile(statepath,stateid,procid):
    resultsfile = _open_file_reading( os.path.join(statepath,str(stateid),"procdata",str("results_%d.dat" % procid) ) )
    lines = resultsfile.readlines()
    resultsfile.close()
    return lines

######################################################################################################

def _open_file_reading(filename):
    fileobject = None
    if(os.path.isfile(filename)):
        try:
            fileobject = open(filename,"r")
        except IOError:
            print "Error: could not read file %s" % filename
    return fileobject

######################################################################################################
        
def _open_file_writing(filename):
    mode = "o"
    if(os.path.isfile(filename)):    
        print "File %s already exists. (o)verwrite/(a)ppend/(q)uit ?" % filename,
        inline = sys.stdin.readline()
        mode = inline.split()[0]
    
    if mode=="o":       fileobject = open(filename,"w")
    elif mode=="a":     fileobject = open(filename,"a")
    else: fileobject =  None
    return fileobject    
