#!/opt/local/bin/python

import conf_convert
import ase
import sys

def view_con(filename='reactant.con', rep=1):
    print filename[-3:]
    if filename[-3:] == 'con':
        a = conf_convert.con2Asap(filename)[0]
        conf_convert.asap_translate_and_wrap(a, [5,5,5])
    else:
        a = ase.read(filename)
    if filename[-3:] == 'xyz':
        rep_in_plane = 1
        print '\n####################################'
        print 'xyz does not support box repetitions'
        print '####################################\n'

    pos = a.get_positions()
    cell = a.get_cell()
    max_z = max(pos[:,2])
    if cell[2][2] > max_z+10:
        # seems to be surface
        rep = [rep,rep,1]
    else:
        rep = [rep,rep,rep]
    
    p = ase.view(a.repeat(rep))
#    raw_input('seen enough?')


if __name__ == "__main__":
    if 2 < len(sys.argv):
        filename = sys.argv[2]
    else:
        filename = 'reactant.con'
    
    if 1 < len(sys.argv):
        rep = int(sys.argv[1])
    else:
        rep = 1
    view_con(filename , rep)
