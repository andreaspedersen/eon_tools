#!/usr/bin/python
import ase

import eon_read

######################################################################################################

def convert_to_ase(passed):
    '''
    Convert structure to ase atoms object
    '''
    eon_atoms = passed['reactant_config']    
    return ase.Atoms(symbols=eon_atoms.names, positions=eon_atoms.r, cell=eon_atoms.box)


def convert_ase_atom_types(ase_atoms, conversions):
    '''
    Convert atom nrs which are given in a list of pairs 
    [[2,4],[5,1],...] will change 2 to 4 and 5 to 1
    '''
    
    atom_nrs_before = ase_atoms.get_atomic_numbers()
    atom_nrs_after = []
    atom_nr_after = None 
    
    for atom_nr_before in atom_nrs_before:
        atom_nr_after = atom_nr_before
        for conversion in conversions:
            if conversion[0] == atom_nr_before:
                atom_nr_after = conversion[1]
                break
        atom_nrs_after.append(atom_nr_after)
    ase_atoms.set_atomic_numbers(atom_nrs_after)



#s = eon_extract.eon_read.read_statelist()
#a = s[0]
#a_a = eon_ase.convert_to_ase(a)
#eon_ase.convert_ase_atom_types(a_a,[[1,29],[2,40]])