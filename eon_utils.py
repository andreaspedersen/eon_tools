#!/usr/bin/python
import os
from optparse import OptionParser
import sys
import atom
from atom import pbc
import numpy
import math

def open_file_reading(filename):
    fileobject = None
    if(os.path.isfile(filename)):
        fileobject = open(filename,"r")
    return fileobject
        
def open_file_writing(filename):
    mode = "o"
    if(os.path.isfile(filename)):    
        print "File %s already exists. (o)verwrite/(a)ppend/(q)uit ?" % filename,
        inline = sys.stdin.readline()
        mode = inline.split()[0]
    
    if mode=="o":       fileobject = open(filename,"w")
    elif mode=="a":     fileobject = open(filename,"a")
    else: fileobject =  None
    return fileobject

def determine_num_states(statepath):
    '''
	Determines the total number of states which the AKMC simulaton has visited. 
	This is simply determined from the number of lines in the state_table file.
    '''
    state_table = open_file_reading(os.path.join(statepath,"state_table"))
    num_states = 0
    if state_table != None:
        lines = state_table.readlines()
        state_table.close()
        num_states = len(lines)
    return num_states

def determine_num_procs(statepath,stateid):
    '''
	Returns the number of unique processes leading out of a specified state.
    '''
    path = os.path.join(statepath,str(stateid) )
    if not os.path.isdir(path):
        print "error: %s is not a valid statepath" % path
    processtable = open_file_reading(os.path.join(path,"processtable"))
    numprocs = 0
    if processtable != None:
        lines = processtable.readlines()
        processtable.close()
        numprocs = len(lines)-1
    return numprocs
    
    
    
    
def read_proc_info(statepath,stateid,procid):
    path = os.path.join(statepath,str(stateid) )
    if not os.path.isdir(path):
        print "error: %s is not a valid statepath" % path
    processtable = open_file_reading(os.path.join(path,"processtable"))
    barrier = None
    procdict = { }
    if processtable !=None:
        lines = processtable.readlines()
        processtable.close()
        assert(procid<len(lines))
        barrier =float( lines[procid+1].split()[6]  )
        prefactor = float(lines[procid+1].split()[2] )
        productid = int(lines[procid+1].split()[3] )
        saddle_energy = float(lines[procid+1].split()[1] )
        product_energy = float(lines[procid+1].split()[4] )
        procdict = {'barrier': barrier,
                    'prefactor': prefactor,
                    'productid': productid,
                    'reactantid': stateid,
                    'saddle_energy': saddle_energy,
                    'product_energy': product_energy
                    }
    return procdict
    
def read_state_energy(statepath,stateid):    

    infofile = open_file_reading( os.path.join ( statepath, str(stateid), "info" ) )
    for line in infofile:
        if len(line.split() ) > 2:
            if line.split()[0] == "reactant" and line.split()[1] == "energy":
                energy = float( line.split()[3] )
                break
    infofile.close()
    return energy


def read_state(statepath,stateid):
    
    reactantcon = atom.Atoms(0)
    
    reactantfilename = os.path.join(statepath,str(stateid),"reactant.con")
    assert(os.path.isfile(reactantfilename) )
    reactantcon.read_con(reactantfilename)
    return reactantcon
    
    
    
def read_proc_con(statepath,stateid,procid):
    processpath = os.path.join(statepath,str(stateid),"procdata")
    reactantcon = atom.Atoms(0)
    saddlecon   = atom.Atoms(0)
    productcon  = atom.Atoms(0)
    reactantcon.read_con(os.path.join(processpath,"reactant_"+str(procid)+".con"))
    saddlecon.read_con(os.path.join(processpath,"saddle_"+str(procid)+".con"))
    productcon.read_con(os.path.join(processpath,"product_"+str(procid)+".con"))
    return  reactantcon,saddlecon,productcon 



def read_full_proc_data(statepath,stateid,procid):
    '''
	Read all data from a specified process into a dictionary and returns it
    '''
    procdict= read_proc_info(statepath,stateid,procid)
    
    reactantcon,saddlecon,productcon = read_proc_con(statepath,stateid,procid)
    procdict['id']=procid
    procdict['reactant']=reactantcon
    procdict['product']=productcon
    procdict['saddle']=saddlecon

    return procdict        



def read_full_state_data(statepath,stateid):
    '''
	Reads all data from a specified state into a dictionary and returns it.
    '''
    print "reading state %d" % stateid
    
    statedict = { }
    statedict['id'] = stateid
    statedict['numprocs'] = determine_num_procs(statepath,stateid)
    statedict['reactant'] = read_state(statepath,stateid)
    statedict['energy'] = read_state_energy(statepath,stateid)
    statedict['proclist'] = [ ]
    for procid in range(statedict['numprocs'] ):
        statedict['proclist'].append ( read_full_proc_data(statepath,stateid,procid) )

    return statedict


def read_full_data(statepath):
    '''
	Read data from all states and return it in a list containing dictionaries with information
    '''
    numstates = determine_num_states(statepath)
    statelist = [ ]
    numstatess=0
    for stateid in range(numstates):
        statelist.append( read_full_state_data(statepath,stateid) )
        numstatess+=1
       # if numstatess > 10:
       #     break
            
    return statelist



def read_superbasin(superbasinpath,sbfilename):
        sbfile = open(os.path.join(superbasinpath,sbfilename),"r")

        line = sbfile.readline().split()
        numstates = len(line)
        states = [ ]
        for i in range(numstates):
            states.append( int(line[i]) )
        sbdict = { 'numstates': numstates, 'states':states}
        sbfile.close()    
        return sbdict     



def read_superbasin_data(superbasinpath):
    '''
	Reads superbasin information and returns a list of all superbasins.
	This list contains dictionaries with further information about the superbasins
    '''
    superbasins = []
    #determine the number of current superbasins
    for filename in os.listdir(superbasinpath):
        if filename == "storage":
            continue
        superbasinid=int(filename)
        sbdict = read_superbasin(superbasinpath,filename)
        sbdict['id'] = superbasinid
        sbdict['old'] = False
        superbasins.append(sbdict)

    #determine all old superbasins
    oldsuperbasins = []
    
    if os.path.isdir( os.path.join(superbasinpath,"storage") ):
        for filename in os.listdir( os.path.join( superbasinpath,"storage") ):
            superbasinid = int(filename)
            sbdict = read_superbasin( os.path.join( superbasinpath,"storage"), filename)
            sbdict['id'] = superbasinid
            sbdict['old'] = True
            oldsuperbasins.append(sbdict)
    
    #merge superbasins into list with good id's at proper place.
    maxid = max( [ i['id'] for i in superbasins+oldsuperbasins ] )
    sblist = [None] *( maxid + 1 )
    for sb in superbasins+oldsuperbasins:
        sblist[ int( sb['id'] ) ] = sb
        
    return sblist



def distance_moved(con1,con2,atomid):
    assert(atomid<con1.num_atoms)
    assert(atomid<con2.num_atoms)
    return numpy.linalg.norm(con1.r[atomid] - con2.r[atomid])


    
def calc_com_co(r_c,r_o,box,ibox):

    weightc = 0.428649280166606
    weighto = 1.0-weightc
    r_com = r_c + weighto*pbc(r_o-r_c,box,ibox)

    return r_com    



def co_rotation(rc1,ro1,rc2,ro2,box,ibox):
    v1 = pbc( ro1-rc1, box,ibox)
    v2 = pbc( ro2-rc2, box,ibox)
    angle = math.acos( numpy.dot(v1,v2)  / ( numpy.linalg.norm(v1) * numpy.linalg.norm(v2) ) ) * 180.0/math.pi
    return angle



def h2o_rotation(rha1,rhb1,ro1,rha2,rhb2,ro2,box,ibox):
    '''
    calculates the rotation of a water molecule between two different configurations.
    To do this we first calculate the vector from O to the midpoint of the two H atoms for both configurations.
    The angle between these vectors is then computed. 
    '''
    pass 
    
def calc_co_translate(statelist,outfile=None,ncomax=None):
    n_molecules,n_h2o,n_co = num_molecules( statelist[0]['proclist'][0]['reactant'] )    
    
    if ncomax==None:
        ncomax=n_co
    
    #write header
    if outfile!=None:
        outfile.write("CO center of mass translations\n")
        outfile.write("#St-id\tProc-id\tBarrier\t")
        for i in range(n_co):
            outfile.write("CO_%s\t" % str(i) ) 
        outfile.write("Average\n")
    
    #calculate displacement
    box = statelist[0]['proclist'][0]['reactant'].box
    ibox = numpy.linalg.inv(box)
    
    for state in statelist:
        for proc in state['proclist']:
            total_disp = 0.0
            if outfile!=None:
                outfile.write("%d\t%d\t%f\t" % (state['id'],proc['id'],proc['barrier'] ) ) 
            for ico in range(min(ncomax,n_co)):
                id_c = n_h2o*3+ico
                id_o = n_h2o*3+n_co+ico

                com_disp = pbc ( calc_com_co( proc['reactant'].r[id_c] , proc['reactant'].r[id_o],box,ibox ) - calc_com_co( proc['product'].r[id_c] , proc['product'].r[id_o],box,ibox ), box,ibox )
                proc['co_disp']=com_disp
                com_disp = numpy.linalg.norm(com_disp)

                total_disp += com_disp 
                if outfile!=None:
                    outfile.write("%f\t" %  com_disp)
            if outfile!=None:                    
                outfile.write("%f\n" % ( total_disp/n_co ) )




def calc_co_rotation(statelist,outfile):
    n_molecules,n_h2o,n_co = num_molecules( statelist[0]['proclist'][0]['reactant'] )    
    
    #write header
    outfile.write("#CO rotation angles\n")
    outfile.write("#Barrier\t")
    for i in range(n_co):
        outfile.write("CO_%s\t" % str(i) ) 
    outfile.write("Average\n")    

    box = statelist[0]['proclist'][0]['reactant'].box
    ibox = numpy.linalg.inv(box)

    for state in statelist:
        for proc in state['proclist']:
            total_angle = 0.0
            outfile.write("%f\t" % proc['barrier'] )
            for ico in range(n_co):
                id_c = n_h2o*3+ico
                id_o = n_h2o*3+n_co+ico
                angle = co_rotation(proc['reactant'].r[id_c],
                                    proc['reactant'].r[id_o],
                                    proc['product'].r[id_c],
                                    proc['product'].r[id_o],
                                    box,
                                    ibox
                                    )
                total_angle +=angle
                outfile.write("%f\t" % angle )
            outfile.write("%f\n" % ( angle / n_co) )

    
    
    
def num_molecules(atomconfig):
    n_h2o = 0
    n_co  = 0

    for iatom in range(atomconfig.num_atoms):
        if atomconfig.names[iatom] == 'H':
            n_h2o+=1
        else:
            break
    n_h2o = int ( n_h2o / 2 )
    n_co = int ( (atomconfig.num_atoms - n_h2o*3 ) / 2 )
    
    return n_co+n_h2o, n_h2o,n_co    
