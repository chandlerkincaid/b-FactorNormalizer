This tool was created from the folowing formula:

Bnorm = (B - Bmean)/Bstdev

Taken from the following paper:

Protein Flexibility and Rigidity Predicted From Sequence
Shlessinger & Rost (2005)

To use save the file, change permission to allow running as executable, and run from terminal. Type the file name
b-FactorNorm.py --help for more information.

This program uses the Biopython, Numpy, and Argparse modules
If you do not have these modules they can be installed via "pip install" or "conda install" for Anaconda users.
Please see their respective documentation for details.

You can now also get a raw csv output for the atoms by adding -w to your command.

You can now also normalize per chain with -c and normalize per residue range(experimental) with -r int int
