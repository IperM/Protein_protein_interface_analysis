#!/usr/bin/env python

"""
    Initial setup a structure for Energy evaluation
"""
import argparse
import sys
import os
import numpy as np

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

print(vars(st[0]['A'][42]['O']))
print(vars(st[0]['A'][42]['O'].xtra['vdw']))


def Eelec(distance, ch1, ch2):
    denominator=1-7.7839*np.exp(-0.3153*distance)

    Er=(86.9525/denominator)-8.5225 # This is a constant called Melher-Solajer Dielectric

    Eelc=332.16*((ch1*ch2)/(Er*distance)) #here we get the electreostatics energy
    return(Eelc)

def Evdw(distance, sig1, sig2, eps1, eps2):
    
    cof1=((((sig2)**6)*(sig1)**6)/((distance)**12))
    cof2=((((sig2)**3)*(sig1)**3)/((distance)**6))

    Evdw=(4*np.sqrt(eps1*eps2))*(cof1-cof2)# here we get van der waals energy
    return(Evdw)

def Esolvation(sig, ASA):
    R = sig+float(ASA)
    return(R)
atoms_chainA = []
atoms_chainE = []
interface_distan = 6.0
Total_Eelec = 0
Total_Vdw = 0
Total_atoms = 0
Total_Solvation = 0  
Solvation_Atom1 = 0
Solvation_Atom2 = 0

for atom1 in st.get_atoms():
    for atom2 in st.get_atoms():
        if atom1.get_parent().get_parent() != atom2.get_parent().get_parent():    	    
            distances = np.linalg.norm(atom1.get_coord()-atom2.get_coord())
            if float(distances) <= interface_distan:
                Total_atoms += 2 
                print('Atom 1 with ID: ',atom1.id[0])
                print('Atom 2 with ID: ',atom2.id[0])
                print('The distance between this 2 atoms is: ',distances)
                   #print(vars(atom1))
                   #print(vars(atom1.xtra['vdw']))
                print(atom1.get_parent().get_parent(),atom2.get_parent().get_parent())
                ch1 = atom1.xtra['charge']
                ch2 = atom2.xtra['charge']
                   #print(ch1,ch2,distances)
                print('Electric Energy of this 2 atoms: ',Eelec(distances, ch1, ch2))
                sig1 = vars(atom1.xtra['vdw'])
                sig1 = sig1['sig']
                sig2 = vars(atom2.xtra['vdw'])
                sig2 = sig2['sig']
                eps1 = vars(atom1.xtra['vdw'])
                eps1 = eps1['eps']
                eps2 = vars(atom2.xtra['vdw'])
                eps2 = eps2['eps']
                print('Van der Waals Energy of this 2 atoms: ',Evdw(distances, sig1, sig2, eps1, eps2))
                sigma = vars(atom1.xtra['vdw'])
                sigma = sigma['fsrf']
                try:
                    Asa =  atom1.xtra['EXP_NACCESS']
                except KeyError:
                    Asa = 0.0
                print(sigma,Asa)
                sigma2 = vars(atom2.xtra['vdw'])
                sigma2 = sigma2['fsrf']
                try:
                    Asa2 =  atom2.xtra['EXP_NACCESS']
                except KeyError:
                    Asa2 = 0.0  
                Solvation_Atom1 = Esolvation(sigma,Asa)
                Solvation_Atom2 = Esolvation(sigma2,Asa2)
                print('Solvation Energy Atom1: ',Esolvation(sigma,Asa))
                print('Solvation Energy Atom2: ',Esolvation(sigma2,Asa2))
                print('Total Solvation Energy of this pair of atoms: ',Solvation_Atom1+Solvation_Atom2)
                   
                print("Pair of atom number: ",Total_atoms)
                Total_Eelec += Eelec(distances,ch1,ch2)
                Total_Vdw += Evdw(distances, sig1, sig2, eps1, eps2)
                Total_Solvation += Solvation_Atom1+Solvation_Atom2
                   
                print('---------------------------------------')
print('Total Electric Energy of the interface: ',Total_Eelec)
print('Total Van der Waals Energy: ',Total_Vdw)
print('Solvation Energy of the interface: ',Total_Solvation)



