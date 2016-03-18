''' 
This scripts splits a large posterior_sample.dat file into small files. This is useful in evolving posterior samples in parallel.

Usage:
python split_posterior_samples.py <posterior_sample_file> <tag name for the small files> <directory path where the small files will be created> <No of
 files> 
'''
#(C) Bhooshan Gadre, Anuradha Gupta, 16-03-2016

from itertools import izip_longest
import argparse

# It will create either no_of_files or no_of_files+1 files

# Set up the parsing
parser = argparse.ArgumentParser(description = 'spliting the posterior samples files into small files.')
parser.add_argument("datfile", help = "The big posterior_samples.dat file to read in")
parser.add_argument("tag", help = "tag for small output files")
parser.add_argument("outdir", help = "directory for output file (default: current directory)", default = ".")
parser.add_argument("no_of_files", help = "Number of smaller files", type=int)
args = parser.parse_args()

datfile = args.datfile
tag = args.tag
no_of_files = args.no_of_files

def grouper(n, iterable, fillvalue=None):
    "Collect data into fixed-length chunks or blocks"
    # grouper(3, 'ABCDEFG', 'x') --> ABC DEF Gxx
    args = [iter(iterable)] * n
    return izip_longest(fillvalue=fillvalue, *args)


with open(datfile) as f:
    heads = f.readline()
    totn = len(f.readlines())

with open(datfile) as f:
    n = int(totn/no_of_files)
    for i, g in enumerate(grouper(n, f, fillvalue=''), 0):
        of = '{0}.{1}'.format(tag, i)
        with open(of, 'w') as fout:
	    if i != 0:
	        fout.writelines(heads)
            fout.writelines(g)
    
