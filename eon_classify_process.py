#!/usr/bin/python
import os
from eon_utils import *
import numpy
from optparse import OptionParser
import time

parser = OptionParser()
parser.add_option("-d", "--dir",type="string", dest="directory", help="Path of the directory eon results", metavar="DIR",default=os.getcwd())
parser.add_option("-o", "--output",type="string", dest="output", help="Directory to write results to")
parser.add_option("-i", "--co_id", type="int",dest="co_id",help="ID of the co molecule to track")
parser.add_option("-c", "--cutoff", type="float",dest="barrier_cutoff",help="Maximum barrier height of a process to be taken into account",default=0.1)
parser.add_option("-p", "--include_procs", action="store_true", dest="include_procs",help = "Also analyze product states leading out of visited states", default=False) 
(options, args) = parser.parse_args()


###########################################################################################
###########################################################################################

def calc_neighbours_r(r,cutoff,dist_array):
    '''
    Calculates a list of all id's in coord_array which are within a cubic
    box of side cutoff from r
    dist_array is an Nx4 numpy array containing the id's of the moilecules in the
    zeroth column followed by the coordinates:
    [[id,x,y,z],
     [id,x,y,z],
     ....
     [id,x,y,z]]
     
    Return value: numpy array with id's which match
    '''
    assert(dist_array.ndim == 2)
    assert(dist_array.shape[1] == 4)
    

    for i in [3,2,1]:
        istart=0
        iend=0    

        dist_array = dist_array[ dist_array[:,i].argsort() ] 
        
        index = 0
        start_found=False
        end_found=False
        for j in range(len( dist_array ) ):
            if ( dist_array[j][i] >= - cutoff and start_found==False ):
                istart = j
                start_found = True
            elif ( dist_array[j][i] >= cutoff and end_found ==False):
                iend = j
                end_found = True
                break
        if start_found==True and end_found==False:
            iend = len(dist_array) 
        if start_found==False:
            print "NO RESULTS"
            return []
        dist_array = dist_array[istart:iend]    
    return numpy.sort( dist_array[:,0] )

###########################################################################################
###########################################################################################


def atomid_to_molid(n_molecules,n_h2o,n_co,atomid):
    #hydrogen in water molecule
    assert(type(atomid) is int)
    if atomid < 2 * n_h2o:
        return atomid/2
   #oxygen in water molecule
    elif atomid < 3 * n_h2o:
        return atomid - 2*n_h2o
    #carbon in CO        
    elif atomid < 3 * n_h2o + n_co:          
        return atomid - 3 * n_h2o
    #oxygen in CO
    elif atomid < 3*n_h2o + 2*n_co:
        return atomid - 3*n_h2o - n_co
    else:
        print "ERROR %d is not a valid atomid" % atomid 

###########################################################################################
###########################################################################################

        
def molid_to_atomids(n_molecules,n_h2o,n_co,molid):
    assert(type(molid) is int)
    if molid < n_h2o:
        return [ molid*2, molid*2+1, molid+2*n_h2o ]
    elif molid <= n_molecules:
        return [ molid + 3*n_h2o, molid + 3*n_h2o + n_co ]

###########################################################################################
###########################################################################################

    
def find_water_near_r(r,cutoff,config):
    '''
    Calculates all oxygen neighbors with specified cutoff of r using sweep and prune
    Returns a sorted list of atomic ids of all atoms in the water molecules found within the cube
    '''
    box = config.box
    ibox = numpy.linalg.inv(box)

    n_molecules, n_h2o, n_co = num_molecules( config )
    
    #numpy array with all oxygen ids.
    oxygenids = numpy.zeros( (n_h2o,1) )
    for i in range(n_h2o): oxygenids[i] = i+n_h2o*2

    #calculate the pbc corrected distance from r to every oxygen atom in the system.
    pbcdist = pbc( config.r - r, box,ibox)    
    
    #numpy array with all oxygen id's and the corresponding distances. need as input for sweep and prune function
    oxygenlist = numpy.concatenate( (oxygenids,pbcdist[2*n_h2o:3*n_h2o]), axis=1)
    oxygenidlist = calc_neighbours_r(r,cutoff,oxygenlist)
    
    #add hydrogen belonging to the oxygens
    neighbourlist = []

    for oxygenid in oxygenidlist:
        molid = atomid_to_molid(n_molecules,n_h2o,n_co,int (oxygenid) )
        neighbourlist += molid_to_atomids(n_molecules,n_h2o,n_co,molid)
    
    
    #now calculate the actual distances and store in a nice array
    nearest_neighbours = numpy.zeros( ( len(neighbourlist), 2 ) )
    for i in range(len(neighbourlist)):
        nearest_neighbours[i][0] = neighbourlist[i]
        nearest_neighbours[i][1] = numpy.linalg.norm(pbcdist[ neighbourlist[i] ] )
    nearest_neighbours = nearest_neighbours [ nearest_neighbours[:,1].argsort() ]
    

    return nearest_neighbours    

###########################################################################################
###########################################################################################
 
def surface_site_type(neighbour_list,config,box,ibox,n_h2o,n_co,r_com):
    
    mindist = box[0][0]    
    numcheck =0 
    molidlist=[]
    distlist = []
    distdict = { }

    neighbour_coordinates = numpy.zeros ( ( len(neighbour_list), 4) )

    for i in neighbour_list:

        if i[0] >= 2*n_h2o and i[0] < 3*n_h2o:
            pbcdist = pbc(config.r[ int( i[0] ) ] - r_com, box,ibox)
            distance =  math.sqrt(pbcdist[0]*pbcdist[0] + pbcdist[1]*pbcdist[1])
            molid = atomid_to_molid( n_h2o+n_co, n_h2o, n_co,int( i[0] ) )
            distdict[distance]= molid            

            neighbour_coordinates[numcheck][0] = i[0]
            for j in range(1,4): neighbour_coordinates[numcheck][j] = config.r[ int(i[0]) ][j-1] 

            mindist = min( mindist, distance )    
            numcheck +=1


        if numcheck >= 6:
            break
    neighbour_coordinates = neighbour_coordinates[0:numcheck]
    
    items = distdict.items()
    items.sort()
    molidlist =  [value for key, value in items]
     

    
    if mindist < 1.0:
        site_hydrogen_type = hydrogen_type(molidlist[1:4],config,box,ibox,n_h2o,n_co,r_com)
        return "crystal", site_hydrogen_type
    elif mindist > 1.7:
        #Sort oxygen atoms by z coordinate
        neighbour_coordinates = neighbour_coordinates [ neighbour_coordinates[:,3].argsort() ]
        molidlist = [ atomid_to_molid(n_h2o+n_co, n_h2o,n_co,int(i) ) for i in neighbour_coordinates[:,0] ]
        site_hydrogen_type = hydrogen_type(molidlist[::-1][0:3], config,box,ibox,n_h2o,n_co,r_com)
        return "chute", site_hydrogen_type
    else:
        return "unknown", "unknown"

###########################################################################################
###########################################################################################

def site_classify_numprotons(num_hydrogens_up):
    if num_hydrogens_up==0:
        return "C site"
    elif num_hydrogens_up==1:
        return "A site"
    elif num_hydrogens_up==2:
        return "B site"
    elif num_hydrogens_up==3:
        return "D site"    
    else:
        return "unknown" 


###########################################################################################
###########################################################################################

def hydrogen_type(molidlist,config,box,ibox,n_h2o,n_co,r_com):

    num_hydrogens_up = 0
    for molid in molidlist:
        atomlist = molid_to_atomids(n_h2o+n_co,n_h2o,n_co,molid)
        assert(len(atomlist)==3)
        #loop over hydrogens
        print atomlist
        for atomid in atomlist[:2]:
            if config.r[atomid][2] > ( config.r[ atomlist[2] ][2] + 0.4 ):
                num_hydrogens_up+=1
    
    return site_classify_numprotons(num_hydrogens_up)

        
###########################################################################################
###########################################################################################


def analyze_config(config,cutoff,co_id):
        
        n_molecules, n_h2o, n_co = num_molecules( config )        
        id_c = n_h2o*3 + ( co_id - n_h2o ) 
        id_o = id_c + n_co
        box = config.box
        ibox = numpy.linalg.inv(box)
        

         #calculate co position
        r_com = calc_com_co( config.r[id_c], config.r[id_o ] , box, ibox )
        neighbour_list = find_water_near_r(r_com,4.2,config)
        site_type, hydrogen_site_type = surface_site_type(neighbour_list,config,box,ibox,n_h2o,n_co,r_com)
        
        return site_type, hydrogen_site_type

###########################################################################################
###########################################################################################
        
        
def analyze_state(state , cutoff, co_id,include_proc):

        #determine site_type
        reactant_site_type, reactant_hydrogen_site_type = analyze_config(state['reactant'],cutoff,co_id)
        state['site_type'] = reactant_site_type
        state['hydrogen_site_type'] = reactant_hydrogen_site_type
        print "\tSite type: %s \n\tHydrogen_site_type: %s\n\n" % (state['site_type'],state['hydrogen_site_type'] )
        
        if include_proc==True:    
            for proc in state['proclist']:
                if proc['barrier'] < 0.1:
                    print "analyzing state %d process %d" % ( state['id'] , proc['id'] ) 
                    site_type,hydrogen_site_type = analyze_config(proc['product'],cutoff,co_id)
                    proc['product_site_type'] = site_type
                    proc['reactant_site_type'] = reactant_site_type
                    proc['product_hydrogen_site_type'] = hydrogen_site_type
                    proc['reactant_hydrogen_site_type'] = reactant_hydrogen_site_type
                    
                    print "\tSite type: %s \n\tHydrogen_site_type: %s\n\n" % (proc['product_site_type'],proc['product_hydrogen_site_type'] )
                    
                    if proc['barrier']<0.0:
                        print " "
                        print "ERROR: state %d, proces %d" % (state['id'],proc['id'] )

###########################################################################################
###########################################################################################


def write_xyz_per_type(statelist,outputdir):
    crystalxyzfile_A = open_file_writing(os.path.join(outputdir,"crystal_sitelist_A.xyz"))
    crystalxyzfile_B = open_file_writing(os.path.join(outputdir,"crystal_sitelist_B.xyz"))
    crystalxyzfile_C = open_file_writing(os.path.join(outputdir,"crystal_sitelist_C.xyz"))
    crystalxyzfile_D = open_file_writing(os.path.join(outputdir,"crystal_sitelist_D.xyz"))            
    crystalxyzfile_unknown = open_file_writing(os.path.join(outputdir,"crystal_sitelist_unknown.xyz"))            


    chutexyzfile_A = open_file_writing(os.path.join(outputdir,"chute_sitelist_A.xyz"))
    chutexyzfile_B = open_file_writing(os.path.join(outputdir,"chute_sitelist_B.xyz"))
    chutexyzfile_C = open_file_writing(os.path.join(outputdir,"chute_sitelist_C.xyz"))
    chutexyzfile_D = open_file_writing(os.path.join(outputdir,"chute_sitelist_D.xyz"))            
    chutexyzfile_unknown = open_file_writing(os.path.join(outputdir,"chute_sitelist_unknown.xyz"))
       
    unknownxyzfile = open_file_writing(os.path.join(outputdir,"unknown_sitelist.xyz"))
        
    for state in statelist:
        if 'site_type' in state.keys():
            if state['site_type'] == "crystal": 
                if state['hydrogen_site_type'] == "A site":
                    state['reactant'].write_xyz(crystalxyzfile_A)                
                if state['hydrogen_site_type'] == "B site":
                    state['reactant'].write_xyz(crystalxyzfile_B) 
                if state['hydrogen_site_type'] == "C site":
                    state['reactant'].write_xyz(crystalxyzfile_C) 
                if state['hydrogen_site_type'] == "D site":
                    state['reactant'].write_xyz(crystalxyzfile_D)                                                         
                if state['hydrogen_site_type'] == "unknown":
                    state['reactant'].write_xyz(crystalxyzfile_unknown)
                     
            if state['site_type'] == "chute":
                if state['hydrogen_site_type'] == "A site":
                    state['reactant'].write_xyz(chutexyzfile_A)                
                if state['hydrogen_site_type'] == "B site":
                    state['reactant'].write_xyz(chutexyzfile_B) 
                if state['hydrogen_site_type'] == "C site":
                    state['reactant'].write_xyz(chutexyzfile_C) 
                if state['hydrogen_site_type'] == "D site":
                    state['reactant'].write_xyz(chutexyzfile_D)                                                         
                if state['hydrogen_site_type'] == "unknown":
                    state['reactant'].write_xyz(chutexyzfile_unknown)   
                                      
            if state['site_type'] == "unknown": state['reactant'].write_xyz(unknownxyzfile)        
            
    crystalxyzfile_A.close()
    crystalxyzfile_B.close()
    crystalxyzfile_C.close()
    crystalxyzfile_D.close()
    crystalxyzfile_unknown.close()
    
    chutexyzfile_A.close()
    chutexyzfile_B.close()
    chutexyzfile_C.close()
    chutexyzfile_D.close()
    chutexyzfile_unknown.close()    
    
                   
    unknownxyzfile.close()

###########################################################################################
###########################################################################################


def write_proc_types(statelist,outputdir):
#    crystal_to_crystal_file = open_file_writing(os.path.join(outputdir,"crystal_to_crystal.dat"))
#    crystal_to_chute_file = open_file_writing(os.path.join(outputdir,"crystal_to_chute.dat"))
#    chute_to_crystal_file = open_file_writing(os.path.join(outputdir,"chute_to_crystal.dat"))
#    chute_to_chute_file = open_file_writing(os.path.join(outputdir,"chute_to_chute.dat"))
                
#    for state in statelist:
#        for proc in state['proclist']:
#            if 'product_site_type' in proc.keys():
#                if proc['product_site_type']=="crystal":
#                    if proc['reactant_site_type']=="crystal":
#                        crystal_to_crystal_file.write("%d\t%d\t%f\n" % (state['id'] , proc['id'] , proc['barrier'] ) )
#                    elif proc['reactant_site_type']=="chute":
#                        crystal_to_chute_file.write("%d\t%d\t%f\n" % (state['id'] , proc['id'] , proc['barrier'] ) )
#                if proc['product_site_type']=="chute":                        
#                    if proc['reactant_site_type']=="crystal":
#                        chute_to_crystal_file.write("%d\t%d\t%f\n" % (state['id'] , proc['id'] , proc['barrier'] ) )
#                    elif proc['reactant_site_type']=="chute":
#                        chute_to_chute_file.write("%d\t%d\t%f\n" % (state['id'] , proc['id'] , proc['barrier']))                        
    
    
#    crystal_to_crystal_file.close()
#    crystal_to_chute_file.close()
#    chute_to_crystal_file.close()
#    crystal_to_chute_file.close()


    proc_classification_file= open_file_writing(os.path.join(outputdir,"proclist_classified.dat") )
    proc_classification_file.write("#reactant\tprocid\t\tproduct\t\tbarrier\t\treac_site\treac_type\tprod_site\tprod_type\n")
    for state in statelist:
        for proc in state['proclist']:
            if 'product_site_type' in proc.keys():
                proc_classification_file.write("%d\t\t%d\t\t%d\t\t%f\t%s\t\t%s\t\t%s\t\t%s\n" % (
                    state['id'],
                    proc['id'],
                    proc['productid'],
                    proc['barrier'],
                    proc['reactant_site_type'].split()[0],
                    proc['reactant_hydrogen_site_type'].split()[0],
                    proc['product_site_type'].split()[0],
                    proc['product_hydrogen_site_type'].split()[0]
                ) )

    proc_classification_file.close()            

###########################################################################################
###########################################################################################
    
                
def main():
    
    #INPUT CHECKS
    
    if "states" not in os.listdir(options.directory):
        print "ERROR: %s is not a valid EON path" % options.directory
        return
    if not options.co_id:
        print "ERROR: must specify id of the co molecule"
        return
    else:
        co_id = int(options.co_id)
    if not options.output:
        print "ERROR: must specify output directory"
        return
    else:
        outputdir=os.path.join(os.getcwd(),options.output)
        if not os.path.isdir(outputdir):
            print "ERROR %s is not a valid directory" % outputdir
            return

    #READ ALL STATES AND PROCESSES
    statepath = os.path.join(options.directory,"states")
    statelist = read_full_data(statepath)

    crystalsitelist = [] 
    chutesitelist   = []
    unknownsitelist = []

    for state in statelist:
        print "Analyzing state %d" % state['id']
        analyze_state(state,4.2, co_id,options.include_procs)


                
    write_xyz_per_type(statelist,outputdir)
    if options.include_procs==True:
        write_proc_types(statelist,outputdir)    
        


###########################################################################################
###########################################################################################

if __name__ == "__main__":
    main()

