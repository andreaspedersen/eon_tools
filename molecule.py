import atom
import re
import sys
import eon_extract
import numpy

class Molecule():
    def __init__(self,name,numatoms,molid,atomids,free=True):
        self._name = name
        self._numatoms = numatoms
        self._molid = molid
        self._atomids = atomids[:]
        self._free = free
        self._com = None
        self._orientation = None

    def __str__(self):
        return "molecule %d (%s)\n\tAtom ids: %s\n\tCOM: %s\n\tOrientation: %s" % (self._molid, self._name, repr(self._atomids), repr(self._com), repr(self._orientation) )


class H2Omol(Molecule):
    #atomid order [H,H,O]
    def __init__(self,molid,atomids,free=True):
        Molecule.__init__(self,"H2O",3,molid,atomids,free)


    def calc_orientation(self,config):
        r1 = eon_extract._fit_box(config.r[ self._atomids[0] ] - config.r[ self._atomids[2] ], config.box)
        r2 = eon_extract._fit_box(config.r[ self._atomids[1] ] - config.r[ self._atomids[2] ], config.box)
        r  = 0.5*(r1+r2)
        r  = r / numpy.linalg.norm(r)
        v = numpy.cross(r1,r2)
        v = v / numpy.linalg.norm(v)
        self._orientation = r           
        self._orientation2 = v     

class COmol(Molecule):
    #atomid order [C,O]
    def __init__(self,molid,atomids,free=True):
        Molecule.__init__(self,"CO",2,molid,atomids,free)

    
    def calc_orientation(self,config):
        r = eon_extract._fit_box(config.r[ self._atomids[1] ] - config.r[self._atomids[0] ], config.box)
        r = r / numpy.linalg.norm(r)
        self._orientation = r

class CO2mol(Molecule):
    #atomid order [C,O,O]
    def __init__(self,molid,atomids,free=True):
        Molecule.__init__(self,"CO2",3,molid,atomids,free)


    def calc_orientation(self,config):
        r1 = eon_extract._fit_box(config.r[ self._atomids[1] ] - config.r[ self._atomids[0] ], config.box)
        r2 = eon_extract._fit_box(config.r[ self._atomids[2] ] - config.r[ self._atomids[0] ], config.box)
        r  = 0.5*(r1-r2)
        r  = r / numpy.linalg.norm(r)
        self._orientation = r


class MolList():
    def __init__(self,config):
        nMolDict = determine_num_molecules(config)
        self._n_h2o = nMolDict['n_h2o']
        self._n_co  = nMolDict['n_co']
        self._n_co2 = nMolDict['n_co2']
        self._num_molecules = self._n_h2o+self._n_co + self._n_co2    

        self._h2oList = []
        self._coList  = []
        self._co2List = []        

        for molid in range(self._n_h2o):
            atomids = [ 2*molid, 2*molid+1, 2*self._n_h2o+molid ] # [ H, H, O]
            free = config.free[atomids[0]] & config.free[atomids[1]] & config.free[atomids[2]]
            mol = H2Omol(molid, atomids,free)
            self._h2oList.append(mol)        
    
        for molid in range(self._n_h2o,self._n_h2o+self._n_co):
            imol  = molid - self._n_h2o
            startid = 3*self._n_h2o
            atomids = [ startid + imol, startid + imol + self._n_co ] # [C,O]
            free = config.free[atomids[0]] & config.free[atomids[1]]
            mol = COmol(molid,atomids,free)
            self._coList.append(mol)

        for molid in range(self._n_h2o+self._n_co,self._n_h2o+self._n_co+self._n_co2):
            imol = molid - self._n_h2o - self._n_co
            startid = 3*self._n_h2o + 2*self._n_co
            atomids = [startid+imol , startid + self._n_co2 + 2*imol, startid + self._n_co2 + 2*imol+1 ] # [C,O,O]
            free = config.free[atomids[0]] & config.free[atomids[1]] & config.free[atomids[2]]
            mol = CO2mol(molid,atomids,free)
            self._co2List.append(mol)
        
        self._list = self._h2oList + self._coList + self._co2List

        #sys.stderr.write("Initialized molcule list: ")
        #sys.stderr.write("H2O: %d, CO: %d, CO2: %d\n" % (self._n_h2o, self._n_co, self._n_co2) )
    
    def __iter__(self):
        for mol in self._list:
            yield mol
    
    def __repr__(self):
        return str(self)

    def __str__(self):
        return "Molecule List contains:\n\tH2O: %d, CO: %d, CO2: %d" % (self._n_h2o, self._n_co, self._n_co2) 

    def __getitem__(self, molid):
        if molid < len(self._list):
            assert (self._list[molid]._molid == molid)
            return self._list[molid]
        return KeyError
        

    def calc_orientations(self,config):
        for mol in self._list:
            mol.calc_orientation(config)

    def h2oList(self):
        for mol in self._h2oList:
           yield mol
           
    def coList(self):
         for mol in self._coList:
            yield mol

    def co2List(self):
        for mol in self._co2List:
            yield mol

    def allMolecule(self):
        for mol in self._list:
            yield mol

    def species(self):
        return ['H2O','CO','CO2']

def determine_num_molecules(config):
    returndict = { }
           #determine the number of molecules
    n_h2o = 0
    n_co  = 0
    n_co2 = 0

    returndict['n_h2o'] = n_h2o
    returndict['n_co']  = n_co
    returndict['n_co2'] = n_co2

    #determine number of h2o's       
    i = 0.0
    for name in config.names:
        if not re.search('^H',name): break
        i+=1
    n_h2o= int(i/2)
    returndict['n_h2o'] = n_h2o
    if len(config)==3*n_h2o:
         return returndict
    ic=0
    io=0
    for name in config.names[3*n_h2o:]:
        if not re.search('^C',name):  break
        ic+=1
    for name in config.names[3*n_h2o+ic:]:
        if ( not re.search('^O',name) ) and ( not re.search('^F',name) ) : break
        io+=1

    if io==ic:
        n_co=io
        returndict['n_co'] = n_co
    elif io==2*ic:
        n_co2 = ic
        returndict['n_co2']= n_co2
        return returndict
    else:
        print ic,io
        print "ERROR: ATOM ORDER IS NOT CORRECT IN REACTANT.CON"
        exit(1)
            
    if len(config)==(3*n_h2o+2*n_co):
        return returndict

    ic =0
    for name in config.names[3*n_h2o+2*n_co]:
        if not re.search('^C',name):  break
        ic+=1
    n_co2 = ic
    returndict['n_co2'] = n_co2
    return returndict


def atomid_to_molid(atomid,molecule_list):
    for molecule in molecule_list:
        for iatom in molecule._atomids:
            if iatom == atomid:
                return molecule._molid
