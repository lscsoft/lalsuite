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

def extract_data(name_of_the_file, arg1, arg2, arg3):
	print  'creating file ' + name_of_the_file + ' (' + arg1 + ' versus ' + arg2 + ' )'
	command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c '+arg1 +' '+arg2 +' '+arg3+ '  > ' + name_of_the_file
	os.system(command)

def extract_data_sum(name_of_the_file, arg1, arg2, arg3):
	print  'creating file ' + name_of_the_file + '( ' + arg1 + ' versus ' + arg2 + ' + ' +arg3+' )'
	command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c '+arg1 +' '+arg2 +' '+arg3+ ' | awk \'{print $2+$3 \" \"$1}\' - > ' + name_of_the_file
	os.system(command)

def extract_data_minus(name_of_the_file, arg1, arg2, arg3):
	print  'creating file ' + name_of_the_file + '( ' + arg1 + ' versus ' + arg2 + ' - ' + arg3+' )'
	command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c '+arg1 +' '+arg2 +' '+arg3+ ' | awk \'{print $2-$3 \" \"$1}\' - > ' + name_of_the_file
	os.system(command)

def Extract2ColAndAdd(name_of_the_file, arg1, arg2):
	print  'creating file ' + name_of_the_file + '( ' + arg1 + ' + ' + arg2 + ') '
	command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c '+arg1 +' ' +arg2 +' | awk \'{print $1+$2}\' - > ' + name_of_the_file
	os.system(command)

def Extract2ColAndSubAnd1Col(name_of_the_file, arg1, arg2, arg3):
    	print  'creating file ' + name_of_the_file + '( ' + arg1 + ' - ' + arg2 +  ' and ' + arg3 + ' )'
	command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c '+arg1 +' ' +arg2 + ' ' + arg3 +' | awk \'{print $1-$2 \" \" $3}\' - > ' + name_of_the_file
	os.system(command)


#starting application here
checkargs()

filename = sys.argv[argc-1]



# first related to the overlap
extract_data_sum('MOverlap.dat', 'overlap' ,'mass1I', 'mass2I')
extract_data('Overlap.dat','overlap','','')
extract_data('OverlapAlphaF.dat','overlap','alpha_f','')
extract_data('OverlapPhase.dat','overlap','phase','')
extract_data('OverlapAlpha.dat','overlap','alpha','')
extract_data('OverlapPsi0T.dat','overlap','psi0T','')
extract_data('OverlapPsi0I.dat','overlap','psi0I','')
extract_data('OverlapPsi3T.dat','overlap','psi3T','')
extract_data('OverlapPsi3I.dat','overlap','psi3I','')
extract_data('OverlapFreqT.dat','overlap','fT','')
extract_data('OverlapFreqI.dat','overlap','fI','')
Extract2ColAndSubAnd1Col('Overlap_deltaF.dat', 'fT', 'fI', 'overlap' )

#related to alpha
extract_data('Alpha.dat','alpha','','')
extract_data('AlphaF.dat','Alpha_f','','')
extract_data('AlphaFpsi0T.dat','Alpha_f','psi0T','')
extract_data('AlphaFpsi0I.dat','Alpha_f','psi0I','')
extract_data('AlphaFpsi3T.dat','Alpha_f','psi3T','')
extract_data('AlphaFpsi3I.dat','Alpha_f','psi3I','')
extract_data('PhaseAlphaF.dat','phase','alpha_f','')
extract_data('PhaseAlpha.dat','phase','alpha','')
extract_data('FreqTAlpha.dat','fT','alpha','')
extract_data('FreqTAlphaF.dat','fT','alpha_f','')
#  absolute value of alpha
command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c overlap alpha_f | awk \'{print sqrt($2*$2) \" \"$1}\' - > OverlapAlphaFAbs.dat'
os.system(command)
extract_data_sum('AlphaFM.dat', 'alpha_f' ,'totalMassT', '')
command = 'lwtprint ' +filename+ ' -t bankefficiency -d \" \" -r 1- -c totalMassT alpha_f | awk \'{print sqrt($1*$1)\" \" sqrt($2*$2)}\' - > AlphaFMAbs.dat'
os.system(command)

#bank related 
extract_data('Psi0IPsi0T.dat','psi0I','psi0T','')
extract_data('Psi3IPsi3T.dat','psi3I','psi3T','')
extract_data('Psi0TPsi3T.dat','psi0T','psi3T','')
extract_data('Psi0IPsi3I.dat','psi0I','psi3I','')


extract_data('OMA.dat','overlap','totalMAssT','alpha_F')

Extract2ColAndAdd('totalMassI.dat', 'mass1I','mass2I')

