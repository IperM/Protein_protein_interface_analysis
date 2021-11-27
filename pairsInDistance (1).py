#!/usr/bin/env python

from Bio.PDB import *
import argparse
import numpy
import sys
import os

parser = argparse.ArgumentParser(
        prog = "pairsInDistance",
        description = "Returns pair residues witht the CA closer than the given distance"
        )
parser.add_argument("file_name", help="PDB file to work with")
parser.add_argument("--distance", help="Max distance",type= int, dest="distance", default= 8)

args = parser.parse_args()

# The file_name and the distance are manddatory so writting them withou "--".
# Storing arguments in variables for a better work:

dis = args.distance
PDBfile = args.file_name

PDBparser = PDBParser(PERMISSIVE = 1) #parse even if some atoms are missing

id_, ext = os.path.splitext(PDBfile)
PDBid = os.path.basename(id_)

structure = PDBparser.get_structure(PDBid, PDBfile)
resListA = structure[0]['A'].get_list()
resListE = structure[0]['E'].get_list()

print("The interface is made upon the pairs of residues distint 6 Ångströms:\n")
for i in resListA:
    for j in resListE:
        
        try:
            D = i['CA'].get_coord() - j['CA'].get_coord()
            distance = numpy.sqrt(numpy.sum(D*D))   #numpy needed as coordinates are arrays
        except KeyError: #There's no CA atom
            continue
        
        if distance <= dis:
            print(i.get_resname(), i.id[1]," in chain A and ", j.get_resname(), j.id[1]," in chain E\n")
                   

