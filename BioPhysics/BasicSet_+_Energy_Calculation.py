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

resListA = st[0]['A'].get_list()
resListE = st[0]['E'].get_list()

Interface_ResiduesA = []
Interface_ResiduesE =[]

for i in resListA:
    for j in resListE:
        
        try:
            D = i['CA'].get_coord() - j['CA'].get_coord()
            distance = np.sqrt(np.sum(D*D))   #numpy needed as coordinates are arrays
        except KeyError: #There's no CA atom
            continue
        
        if distance <= 20:

            if not i.id[1] in Interface_ResiduesA:
                Interface_ResiduesA.append(i.id[1])

            if not j.id[1] in Interface_ResiduesE:
                Interface_ResiduesE.append(j.id[1]) 

            

print("in Chain A ", Interface_ResiduesA)
print("in chain E ", Interface_ResiduesE) 

# First of all we need to define the different variables that we are going to use, they will store the different
# Energy values

interface_distan = 8
Total_Eelec = 0
Total_Vdw = 0
Total_atoms = 0
Total_Solvation = 0  
Solvation_Atom1 = 0
Solvation_Atom2 = 0

print()
print()
print('PAIRS OF ATOMS INSIDE THE INTERFACE AREA OF THE COMPLEX')
print('----------------------------------------------------------------')

# Since we have found the Interface Residues with the script that we have before with 8 Angstroms, we can reduce the search and
# only iterate and compute the different energies from the atoms/resiudes in the Interface, this will make the program run faster 

for i in range(len(Interface_ResiduesA)):
    for j in range(len(Interface_ResiduesE)):  
        for atom1 in st[0]['A'][Interface_ResiduesA[i]].get_atoms():
            for atom2 in st[0]['E'][Interface_ResiduesE[j]].get_atoms():
                if atom1.serial_number < atom2.serial_number: # to check that the atom serial number is always smaller in the first atom than the second 
                    distances = np.linalg.norm(atom1.get_coord()-atom2.get_coord()) # get the distance of the atoms 
                    if distances <= interface_distan: # compare it with the default value that it is 8 Angstroms
                        Total_atoms += 2 # in order to compute 
                        print('Atom 1 with ID: ',atom1.id[0]) # print Atom 1 id
                        print('Atom 2 with ID: ',atom2.id[0]) # print Atom 2 id 
                        print('The distance between this 2 atoms is: ',distances) # print distance between the 2 atoms
                    
                        chain1 = atom1.get_parent().get_parent() # we store the chain of the atom 1 in a variable called chain1
                        chain2 = atom2.get_parent().get_parent() # we store the chain of the atom 2 in a variable called chain2
                        print('Atom 1 chain: ',chain1.id) # print id of chain1 
                        print('Atom 2 chain: ',chain2.id) # print id of chain2
                        ch1 = atom1.xtra['charge'] # we get the atom1 charge that is store inside .xtra dictionary computed before
                        ch2 = atom2.xtra['charge'] # we get the atom2 charge that is store inside .xtra dictionary computed before
                        print('Charge of Atom 1: ',ch1) # print atom1 charge
                        print('Charge of Atom 2: ',ch2) # print atom2 charge
                        print('Electric Energy of this 2 atoms: ',Eelec(distances, ch1, ch2)) # print the Electric Energy of a pair of atoms (atom1 and atom2)
                        sig1 = atom1.xtra['vdw'].sig # as before inside the .xtra dictionary we can get all Van der Wals energies and take the sigma value for atom1
                        sig2 = atom2.xtra['vdw'].sig # as before inside the .xtra dictionary we can get all Van der Wals energies and take the sigma value for atom2
                        eps1 = atom1.xtra['vdw'].eps # as before inside the .xtra dictionary we can get all Van der Wals energies and take the epsilon value for atom1
                        eps2 = atom2.xtra['vdw'].eps # as before inside the .xtra dictionary we can get all Van der Wals energies and take the sigma value for atom2
                    
                        print('Van der Waals Energy of this 2 atoms: ',Evdw(distances, sig1, sig2, eps1, eps2)) # print the Van der Wals Energy of a pair of atoms (atom1 and atom2)
                        sigma = atom1.xtra['vdw'].fsrf # now in order to compute the Solvation Energy of the complex we need the surface value stored in the .xtra dictionary inside fsrf of atom 1
                        if 'EXP_NACCESS' in atom1.xtra: # but some atoms as Hydrogen they don't have fsrf in .xtra dictionary so to avoid a KeyError we need to put the Asa1 variable equal 0 directly
                            Asa = atom1.xtra['EXP_NACCESS']
                        else:
                            Asa = 0.0 # For atoms with no surface value (no .xtra['EXP_NACCESS']
                
                        sigma2 = atom2.xtra['vdw'].fsrf # now in order to compute the Solvation Energy of the complex we need the surface value stored in the .xtra dictionary inside fsrf of atom 2
                        if 'EXP_NACCESS' in atom2.xtra:
                            Asa2 = atom2.xtra['EXP_NACCESS'] # but some atoms as Hydrogen they don't have surface value in .xtra dictionary so to avoid a KeyError we need to put the Asa2 variable equal 0 directly
                        else:
                            Asa2 = 0.0 # For atoms with no surface value (no .xtra['EXP_NACCESS'])
                    
                    
                        Solvation_Atom1 = Esolvation(sigma,Asa) # store the value returned from Esolvation function in a variable called Solvation_atom1
                        Solvation_Atom2 = Esolvation(sigma2,Asa2) # store the value returned from Esolvation function in a variable called Solvation_atom2
                        print('Solvation Energy Atom1: ',Esolvation(sigma,Asa))
                        print('Solvation Energy Atom2: ',Esolvation(sigma2,Asa2))
                        print('Total Solvation Energy of this pair of atoms: ',Solvation_Atom1+Solvation_Atom2)
                        
                        print("Pair of atom number: ",Total_atoms)
                        Total_Eelec += Eelec(distances,ch1,ch2) 
                        Total_Vdw += Evdw(distances, sig1, sig2, eps1, eps2)
                        Total_Solvation += Solvation_Atom1+Solvation_Atom2
                        
                        print('----------------------------------------------------------------')

# now we need to compute bot Solvation Energies from each monomer of the complex (A and E)
Total_SolvationA = 0
for chain in st[0].get_chains():
    print(chain.id)
    if chain.id == 'A':
        for a in range(len(Interface_ResiduesA)):
            for atom1 in st[0]['A'][Interface_ResiduesA[a]].get_atoms():
                print('Atom with ID: ',atom1.id[0])
                chain = atom1.get_parent().get_parent()
                print('Atom chain A')
                sigma = atom1.xtra['vdw'].fsrf # now in order to compute the Solvation Energy of the complex we need the surface value stored in the .xtra dictionary inside fsrf of atom 1
                if 'EXP_NACCESS' in atom1.xtra:
                    Asa = atom1.xtra['EXP_NACCESS'] # for atoms that has the surface computed 
                else:
                    Asa = 0.0 # but some atoms as Hydrogen they don't have surface value in .xtra dictionary so to avoid a KeyError we need to put the Asa2 variable equal 0 

                Solvation_Atom = Esolvation(sigma,Asa) # we store the value returned by the Esolvation function in a variable
                print('Total Solvation Energy of this pair of atoms: ',Solvation_Atom)
            
                Total_SolvationA += Solvation_Atom      
                print('----------------------------------------------------------------')
        print('Total Solvation of Monomer A: ',Total_Solvation) # here we print the total value of the Solvation Energy of the Monomer A

    else:
        Total_SolvationE = 0 # here we reset it to 0 to compute again the Total Solvation of the other monomer
        for n in range(len(Interface_ResiduesE)):
            for atom1 in st[0]['E'][Interface_ResiduesE[n]].get_atoms():
                print('Atom with ID: ',atom1.id[0])
                print('Atom chain E')
                sigma = atom1.xtra['vdw'].fsrf
                if 'EXP_NACCESS' in atom1.xtra:
                    Asa = atom1.xtra['EXP_NACCESS']
                else:
                    Asa = 0.0

                Solvation_Atom = Esolvation(sigma,Asa)
                print('Total Solvation Energy of this pair of atoms: ',Solvation_Atom)
            
                Total_SolvationE += Solvation_Atom      
                print('----------------------------------------------------------------')

            print('Total Solvation of Monomer E: ',Total_SolvationE)

print('Total Electric Energy of the interface: ',Total_Eelec)
print('Total Van der Waals Energy of the interface: ',Total_Vdw)
print('Total Solvation Energy of the interface: ',Total_Solvation)
print('Total Solvation of Monomer A: ',Total_SolvationA)
print("Total ΔG of the complex' interface: ",Total_Eelec + Total_Vdw + Total_Solvation - Total_SolvationA -Total_SolvationE)

# After the program makes everything, it will print the final Energies(Electric, Vdw and Solvation) and it will compute also the 
# ΔG of the complex and print it
