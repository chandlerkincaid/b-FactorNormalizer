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
file_parser.add_argument("-c", "--bychain", help="adding this flag will normalize per chain", action="store_true")
file_parser.add_argument("-r", "--resrange", nargs=2, type=int, help="adding this flag followed by two numbers will "
                                                            "specify only that residue range. Experimental")
args = file_parser.parse_args()
# biopython file handling
PDB_parser = PDBParser()
# structure = PDB_parser.get_structure('myPDB', open(args.pdb_in, "r"))  # get our structure
structure = PDB_parser.get_structure('myPDB', args.pdb_in)  # get our structure
if args.bychain:
    # this option normalizes bFactor per chain instead of the entire structure
    for model in structure:
        for chain in model:
            current_chain = []
            if args.resrange is not None:
                mod_chain = list(chain)[args.resrange[0]: args.resrange[1]]
            else:
                mod_chain = chain
            for residue in mod_chain:
                for atom in residue:
                    current_chain.append(atom.get_bfactor())
            bFactors_mean = np.mean(current_chain)
            bFactors_stDev = np.std(current_chain)
            [x.set_bfactor((x.get_bfactor() - bFactors_mean) / bFactors_stDev) for x in chain.get_atoms()]
else:
    #  this option normalizes bFactor across the entire structure
    if args.resrange is not None:
        mod_res_list = []
        for model in structure:
            for chain in model:
                mod_res_list.append(list(chain)[args.resrange[0]: args.resrange[1]])
        atom_bFactors = [atom.get_bfactor() for residue in mod_res_list for atom in residue]
    else:
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
if args.resrange is not None:
    chain_output = [list(chain)[args.resrange[0]: args.resrange[1]] for chain in structure.get_chains()]
else:
    chain_output = structure.get_chains()
# add an entire chains residues, atoms, and b -factor to rows
pattern = re.compile('[A-Z0-9]+')
for chains in chain_output:
    atom_strings = [re.findall(pattern, str(atoms.get_full_id())) for residues in chains for atoms in residues]
    chain_column = [atoms[2] for atoms in atom_strings]
    residue_column = [atoms[3] for atoms in atom_strings]
    atom_column = [atoms[4] for atoms in atom_strings]
    csv_output.append(chain_column)
    csv_output.append(residue_column)
    csv_output.append(atom_column)
    csv_output.append([atoms.get_bfactor() for residues in chains for atoms in residues])
# transpose matrix to desired format
transposed = [list(row) for row in itertools.zip_longest(*csv_output)]
if args.csvwrite:
    with open(args.output_name + '.csv', 'w') as csvfile:
        mywriter = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        [mywriter.writerow(row) for row in transposed]

