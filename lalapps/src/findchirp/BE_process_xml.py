"""
python script to process xml file given by the BankEfficiency code
Intensive use of awk and lwtprint. 
"""

__author__ = 'Thomas Cokelaer <Thomas.Cokelaer@astro.cf.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'


import sys
import os

argc =  len(sys.argv)


# a help function
def usage():
    print 'you have to provide only one argument which is the name of the file to process'
    sys.exit(1)

# a checking function to parse parameters
def checkargs():
    if argc != 2:
        usage()

#starting application here
checkargs()

filename = sys.argv[argc-1]



print 'creating file MOverlap.dat (Total injected mass versus Overlap)'
command = "lwtprint "+filename+" -t bankefficiency -d \" \" -r 1- -c overlap mass1I mass2I | awk '{print $2+$3 \" \"$1}' - > MOverlap.dat"
os.system(command)

print 'creating file Overlap.dat (Overlap)'
command = "lwtprint "+filename+" -t bankefficiency -d \" \" -r 1- -c overlap  > Overlap.dat"
os.system(command)


print 'creating file Overlap.dat (Overlap versus Alpha_f )'
command = "lwtprint "+filename+" -t bankefficiency -d \" \" -r 1- -c overlap  alpha_f | awk '{print $1 \" \"$1}' - > OverlapAlpha.dat"
os.system(command)



