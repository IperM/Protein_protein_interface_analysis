#!/usr/bin/env python

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os
import numpy as np
from Energyevaluation import *

from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.NACCESS import NACCESS_atomic
from forcefield import VdwParamset

parser = argparse.ArgumentParser(
    prog='structure_setup',
    description='basic structure setup'
)

parser.add_argument(
    '--naccess',
    action='store',
    dest='naccess_bin',
    default=os.path.dirname(os.path.abspath(__file__)) + '/soft/NACCESS/naccess',
    help='Vdw parameters'
)
parser.add_argument(
    '--vdw',
    action='store',
    dest='vdwprm_file',
    default=os.path.dirname(os.path.abspath(__file__)) + '/data/vdwprm',
    help='Vdw parameters'
)

parser.add_argument('pdb_file',help='Input PDB', type=open)
parser.add_argument('pdbqt_file',help='Input PDBQT', type=open)

args = parser.parse_args()

print("PDB File:", args.pdb_file.name)
print("PDBQT File:", args.pdbqt_file.name)
print("Parameters File:", args.vdwprm_file)

# loading VdW parameters
ff_params = VdwParamset(args.vdwprm_file)

parser = PDBParser(PERMISSIVE=1)
print('Parsing PDB', args.pdb_file.name)

# load structure from PDB file of PDB ifle handler
st = parser.get_structure('STR', args.pdb_file.name)

# We will use the xtra attribute in Bio.PDB.Atom to hold the new data
# Getting Charges and Atom type from PDBQT
print('Parsing PDBQT', args.pdbqt_file.name)
params=[{}]

#Fix aton numbers when they do not start in 1
i = 1
for at in st.get_atoms():
    at.serial_number = i
    i += 1

for line in args.pdbqt_file:
    line = line.rstrip()
    #Skip TER records from PDBQT
    if line.find('TER') != -1:
        continue
    params.append({'charge': line[69:76], 'type': line[77:].replace(' ','')})

total_charge = 0.
for at in st.get_atoms():
    at.xtra['atom_type'] = params[at.serial_number]['type']
    at.xtra['charge'] = float(params[at.serial_number]['charge'])
    at.xtra['vdw'] = ff_params.at_types[at.xtra['atom_type']]
    total_charge += at.xtra['charge']
print('Total Charge: {:8.2f}'.format(total_charge))

# Calculating surfaces
# Srf goes to .xtra['EXP_NACCESS'] field
srf = NACCESS_atomic(st[0], naccess_binary=args.naccess_bin)

# Simple Test for atom fields. Choose one atom from the PDB file. st[model][chain_id][res_number][atom_name]

#print(vars(st[0]['A'][42]['O']))
#print(vars(st[0]['A'][42]['O'].xtra['vdw']))



atoms_chainA = []
atoms_chainE = []
interface_distan = 6.0
Total_Eelec = 0
Total_Vdw = 0
Total_atoms = 0
Total_Solvation = 0  
Solvation_Atom1 = 0
Solvation_Atom2 = 0
print()
print('PAIRS OF ATOMS INSIDE THE INTERFACE AREA OF THE COMPLEX')
print('----------------------------------------------------------------')
for atom1 in st.get_atoms():
    for atom2 in st.get_atoms():
        if atom1.get_parent().get_parent() != atom2.get_parent().get_parent():    	    
            distances = np.linalg.norm(atom1.get_coord()-atom2.get_coord())
            if float(distances) <= interface_distan:
                Total_atoms += 2 
                print('Atom 1 with ID: ',atom1.id[0])
                print('Atom 2 with ID: ',atom2.id[0])
                print('The distance between this 2 atoms is: ',distances)
                
                chain1 = atom1.get_parent().get_parent()
                chain2 = atom2.get_parent().get_parent() 
                print('Atom 1 chain: ',chain1.id)
                print('Atom 2 chain: ',chain2.id)
                ch1 = atom1.xtra['charge']
                ch2 = atom2.xtra['charge']
                print('Charge of Atom 1: ',ch1)
                print('Charge of Atom 2: ',ch2)
                print('Electric Energy of this 2 atoms: ',Eelec(distances, ch1, ch2))
                sig1 = atom1.xtra['vdw'].sig
                sig2 = atom2.xtra['vdw'].sig
                eps1 = atom1.xtra['vdw'].eps
                eps2 = atom2.xtra['vdw'].eps
               
                print('Van der Waals Energy of this 2 atoms: ',Evdw(distances, sig1, sig2, eps1, eps2))
                sigma = atom1.xtra['vdw'].fsrf
                if 'EXP_NACCESS' in atom1.xtra:
                    Asa = atom1.xtra['EXP_NACCESS']
                else:
                    Asa = 0.0
         
                sigma2 = atom2.xtra['vdw'].fsrf
                if 'EXP_NACCESS' in atom2.xtra:
                    Asa2 = atom2.xtra['EXP_NACCESS']
                else:
                    Asa2 = 0.0
                
                
                Solvation_Atom1 = Esolvation(sigma,Asa)
                Solvation_Atom2 = Esolvation(sigma2,Asa2)
                print('Solvation Energy Atom1: ',Esolvation(sigma,Asa))
                print('Solvation Energy Atom2: ',Esolvation(sigma2,Asa2))
                print('Total Solvation Energy of this pair of atoms: ',Solvation_Atom1*Solvation_Atom2)
                   
                print("Pair of atom number: ",Total_atoms)
                Total_Eelec += Eelec(distances,ch1,ch2)
                Total_Vdw += Evdw(distances, sig1, sig2, eps1, eps2)
                Total_Solvation += Solvation_Atom1+Solvation_Atom2
                   
                print('----------------------------------------------------------------')
print('Total Electric Energy of the interface: ',Total_Eelec)
print('Total Van der Waals Energy: ',Total_Vdw)
print('Solvation Energy of the interface: ',Total_Solvation)



