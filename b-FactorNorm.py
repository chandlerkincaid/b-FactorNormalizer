#!/usr/bin/env python

import argparse
from Bio.PDB import *
import numpy as np
import os.path

# command line input
file_parser = argparse.ArgumentParser(description="Hello, this tool normalizes pdb b-factors")
file_parser.add_argument("pdb_in", help="specify absolute path and name for input.pdb")
file_parser.add_argument("output_name", help="specify file output name, if no path is included the "
                                             "file will output where the script is located")
args = file_parser.parse_args()
# biopython file handling
PDB_parser = PDBParser()
structure = PDB_parser.get_structure('myPDB', open(args.pdb_in, "r"))  # get our structure
atom_bFactors = [x.get_bfactor() for x in structure.get_atoms()]
bFactors_mean = np.mean(atom_bFactors)
bFactors_stDev = np.std(atom_bFactors)
# set atoms to have normalized bfactor via list comprehensions
[x.set_bfactor(x.get_bfactor() - bFactors_mean / bFactors_stDev) for x in structure.get_atoms()]
# output file
io = PDBIO()
io.set_structure(structure)
io.save(args.output_name)
