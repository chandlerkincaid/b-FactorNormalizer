#!/usr/bin/env python

import argparse
from Bio.PDB import *
import numpy as np
import re
import csv
import itertools

# command line input
file_parser = argparse.ArgumentParser(description="Hello, this tool normalizes pdb b-factors")
file_parser.add_argument("pdb_in", help="specify absolute path and name for input.pdb")
file_parser.add_argument("output_name", help="specify file output name, if no path is included the "
                                             "file will output where the script is located")
file_parser.add_argument("-w", "--csvwrite", help="adding this flag will output a csv file", action="store_true")
args = file_parser.parse_args()
# biopython file handling
PDB_parser = PDBParser()
structure = PDB_parser.get_structure('myPDB', open(args.pdb_in, "r"))  # get our structure
atom_bFactors = [x.get_bfactor() for x in structure.get_atoms()]
bFactors_mean = np.mean(atom_bFactors)
bFactors_stDev = np.std(atom_bFactors)
# set atoms to have normalized bfactor via list comprehensions
[x.set_bfactor((x.get_bfactor() - bFactors_mean) / bFactors_stDev) for x in structure.get_atoms()]
# output file
io = PDBIO()
io.set_structure(structure)
io.save(args.output_name)
# csv output
csv_output = []
chain_output = structure.get_chains()
# add an entire chains residues, atoms, and b -factor to rows
pattern = re.compile('[A-Z0-9]+')
for chains in chain_output:
    atom_strings = [re.findall(pattern, str(atoms.get_full_id())) for atoms in chains.get_atoms()]
    chain_column = [atoms[2] for atoms in atom_strings]
    residue_column = [atoms[3] for atoms in atom_strings]
    atom_column = [atoms[4] for atoms in atom_strings]
    csv_output.append(chain_column)
    csv_output.append(residue_column)
    csv_output.append(atom_column)
    csv_output.append([atoms.get_bfactor() for atoms in chains.get_atoms()])
# transpose matrix to desired format
transposed = transposed = [list(row) for row in itertools.zip_longest(*csv_output)]
if args.csvwrite:
    with open(args.output_name + '.csv', 'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        [mywriter.writerow(row) for row in transposed]

