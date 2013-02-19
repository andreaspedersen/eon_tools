#!/usr/bin/python
import os
from optparse import OptionParser
import sys
import math
import shutil


parser = OptionParser()
parser.add_option("-d", "--dir",type="string", dest="directory", help="Path of the directory eon results", metavar="DIR",default=os.getcwd())
parser.add_option("-o", "--output",type="string", dest="output", help="Directory to write results to")
parser.add_option("-T", "--Temperature", type="float",dest="Tnew",help="New simulation temperature")

(options, args) = parser.parse_args()

###########################################################################################
###########################################################################################


def write_new_processtables(indir,outdir,oldT,newT):
    
    instatesdir = os.path.join(indir,"states")
    outstatesdir = os.path.join(outdir,"states")

    factor = math.exp(oldT/newT)
    processtable_line = "%7d %16.5f %11.5e %9d %16.5f %17.5e %8.5f %12.5e %7d\n"

    for statedir in os.listdir(instatesdir):
        if statedir!="state_table":
            inprocesstablefilename = os.path.join(instatesdir,statedir,"processtable")
            outprocesstablefilename = os.path.join(outstatesdir,statedir,"processtable")
            if os.path.isfile(inprocesstablefilename):
                
                
                inprocesstablefile = open(inprocesstablefilename,"r")
                outprocesstablefile = open(outprocesstablefilename,"w")
                
                #write head line
                line = inprocesstablefile.readline()
                outprocesstablefile.write(line)
                
                for line in inprocesstablefile:              
                    line = line.split()
                    ID = int(line[0])
                    ENERGY = float(line[1])
                    PREFACTOR = float(line[2])
                    PRODUCT = int(line[3])
                    PRODUCT_ENERGY = float(line[4])
                    PRODUCT_PREFACTOR =  float(line[5])
                    BARRIER = float(line[6])
                    RATE = PREFACTOR*math.exp(-BARRIER* 11604.5 / newT )
                    REPEATS = int(line[8])
                    outprocesstablefile.write( processtable_line % (    ID, 
                                                                        ENERGY, 
                                                                        PREFACTOR, 
                                                                        PRODUCT, 
                                                                        PRODUCT_ENERGY, 
                                                                        PRODUCT_PREFACTOR, 
                                                                        BARRIER, 
                                                                        RATE, 
                                                                        REPEATS) )

                inprocesstablefile.close()
                outprocesstablefile.close()
            else:
                print "ERROR: could not open processtable %s" % inprocesstablefilename

###########################################################################################
###########################################################################################


def copyconfigfile(indir,outdir,newT):
    inconfigfilename = os.path.join(indir,"config.ini")
    if os.path.isfile(inconfigfilename):
        inconfigfile = open(inconfigfilename,"r")   
        
    outconfigfilename = os.path.join(outdir,"config.ini")
    outconfigfile = open(outconfigfilename,"w")
    
    for line in inconfigfile:
        if 'temperature' in line.split("=")[0].lower():
            outconfigfile.write("temperature = %f\n" % newT)
        else:
            outconfigfile.write(line)

    inconfigfile.close()
    outconfigfile.close()            

###########################################################################################
###########################################################################################


def get_oldT(indir):
    configfilename = os.path.join(indir,"config.ini")
    if os.path.isfile(configfilename):
        configfile = open(configfilename,"r")
        for line in configfile:
            if 'temperature' in line.split("=")[0].lower():
                oldT = float( line.split("=")[1] )
                break
                
        configfile.close()
    else:
        print "ERROR: could not open config.ini in %s" % indir   
        sys.exit(1)
    
    if not oldT:
        print "ERROR: could not find old temperature in config.ini"
        sys.exit(1)
    
    return oldT

###########################################################################################
###########################################################################################


def main():
    #first check input options
    indir = os.path.join(os.getcwd(),options.directory)
    if not os.path.isdir(indir):
        print "ERROR: input directory %s doesn not exist" % options.directory
        sys.exit(1)
    outdir = os.path.join(os.getcwd(),options.output)
    if not os.path.isdir(outdir):
        print "ERROR: output directory %s doesn not exist" % options.directory
        sys.exit(1)
    if indir==outdir:
        print "ERROR: input directory is the same as the output directory"
        sys.exit(1)
    if len(os.listdir(outdir)) > 0:
        print "ERROR: output directory is not empty"
        sys.exit(1)
    if not options.Tnew:
        print "ERROR: no new temperature specified"
        sys.exit(1)
        
    newT=options.Tnew
    assert(newT > 0.0)

    #Determine old temperature
    oldT = get_oldT(indir)                    
    print "oldT = %f\nnewT = %f\n" %(oldT,newT)
    
    
    #copy reactant.con (copy2() preserves file attributes)
    if os.path.isfile(os.path.join(indir,"reactant.con")):
        shutil.copy2(os.path.join(indir,"reactant.con"),os.path.join(outdir,"reactant.con"))
    
    #copy config.ini, with modified temperature
    copyconfigfile(indir,outdir,newT)            
    
    #copy entire states directory excluding any file named processtable
    instatesdir = os.path.join(indir,"states")
    outstatesdir = os.path.join(outdir,"states")
    shutil.copytree(instatesdir,outstatesdir,symlinks=True,ignore = shutil.ignore_patterns('processtable','superbasin'))
    
    #write new processtables with new rates.
    write_new_processtables(indir,outdir,oldT,newT)
    
    



###########################################################################################
###########################################################################################

if __name__ == "__main__":
    main()

