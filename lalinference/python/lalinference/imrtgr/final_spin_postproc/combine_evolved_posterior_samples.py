'''
This script combines the evolved posterior samples into a big file.

Usage:
python combine_evolved_posterior_samples.py <tag name of the small files> <directory path where the big file will be created> <output file name> <numb
er of files to be combined>
'''
#(C) Anuradha Gupta, 18-03-2016

import argparse
import numpy as np

# Set up the parsing
parser = argparse.ArgumentParser(description = 'combining the posterior sample files into one file.')
parser.add_argument("tag", help = "name tag for input files")
parser.add_argument("outdir", help = "directory for output file (default: current directory)", default = ".")
parser.add_argument("outfilename", help = "output file name")
parser.add_argument("no_of_files", help = "Number of smaller files", type=int)
args = parser.parse_args()

prefix = args.tag
outdir = args.outdir
outfile = args.outfilename
no_of_files = args.no_of_files

# Creating a list of input file names
filenames = [prefix+'_{0}.dat'.format(i) for i in range(no_of_files)]

# Reading in the data from all the files
file_arrays = [np.loadtxt(f) for f in filenames]

# Combining the data 
final_file_array = np.concatenate(file_arrays)

with open(filenames[0]) as f:                                                                   
    h = f.readline()

# saving the data in a new file
np.savetxt(outdir+'/'+outfile, final_file_array, header=h.lstrip('#'))

