"""
python script to create a condor_script and xml prototype  
arguments might be changed directly in this file.
"""

__author__ = 'Thomas Cokelaer <Thomas.Cokelaer@astro.cf.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'


import sys
import os


# init
path  		= '/home/cokelaer/Dev/lalapps/src/findchirp/'
executable_name = 'lalapps_BankEfficiency'
# 

min_mass       	=  3  	   	# miminal individaul mass to inject
max_mass       	=  20.  	# maximal individua lmass to inject
min_massb      	=  3.  		# maximal individua lmass to inject
max_massb      	=  20.  	# maximal individua lmass to inject
signal_psi0_min	= 10000
signal_psi0_max	= 550000
signal_psi3_min	= -4000
signal_psi3_max	= -10
bank_psi0_min	= 10000
bank_psi0_max	= 550000
bank_psi3_min	= -4000
bank_psi3_max	= -10
minimal_match 	=  0.95     	# minimal match in bank creation

lower_frequency = 70		# lower freqency cutoff. identical for both signal and templates for the time being

signal 		= 'EOB'
signal_order 	= 4 
template 	= 'BCV'

#for BCV
fend_min	= -2  
fend_max	= 6 
number_fcut	= 6 
alpha_bank	= 0.01
freq_moment_bank= 2047


noisemodel	= 'LIGOI'
simulation_type = 0       	# 0 for signal only 1 for noise only and 2 for noise+signal
number_of_jobs 	= 100 
signal_amplitude = 0 
number_of_trials_per_job = 20

others = '--print-result-xml  --sampling  4096 --bank-grid-type square'
#others=''


space 		= ' '
# a help function
def usage():
    print 'you have to provide only one argument which is the name of the file to process'
    sys.exit(1)

# a checking function to parse parameters
def checkargs():
    if argc != 2:
        usage()

#starting application here
#checkargs()




arguments = ' --n	' 		+ str(number_of_trials_per_job) \
         +' --mass-range ' 	 + str(min_mass)+ space +str(max_mass) \
         +' --mass-range-bank '	 + str(min_massb)+ space +str(max_massb) \
       	 +' --mm ' 		 + str(minimal_match) \
         +' --fl ' 		 + str(lower_frequency) \
         +' --signal ' 		 + str(signal) \
         +' --signal-order ' 	 + str(signal_order) \
         +' --template ' 	 + str(template) \
         +' --simulation-type '	 + str(simulation_type) \
         +' --signal-amplitude ' + str(signal_amplitude) \
         +' --alpha-bank ' 	 + str(alpha_bank) \
	 +' --noise-model ' 	 + noisemodel \
         +' --fend-bcv ' 	 + str(fend_min) + space + str(fend_max) \
         +' --number-fcut ' 	 + str(number_fcut) \
         +' --bank-psi0-range '	 + str(bank_psi0_min) + space + str(bank_psi0_max) \
         +' --signal-psi0-range '+ str(signal_psi0_min) + space + str(signal_psi0_max) \
         +' --bank-psi3-range '	 + str(bank_psi3_max) + space + str(bank_psi3_min)\
         +' --signal-psi3-range '+ str(signal_psi3_max) + space + str(signal_psi3_min)\
         +' --freq-moment-bank ' + str(freq_moment_bank)\
	 +space + others


# create the condor file
fp =open('BankEfficiency_condor_submit_job','w');
fp.write('Executable   = ' + path + executable_name +'\n')
fp.write('Universe     = vanilla\n')
fp.write('Environment  = LD_LIBRARY_PATH=/software/geopptools/lib\n')
fp.write('Requirements = Memory >=128 && OpSys == "LINUX" && FileSystemDomain == "explorer" && UidDomain == "explorer"\n')
fp.write('+MaxHours =40\n\n')
fp.write('Arguments = '+ '--seed $(Process) ' + arguments)
fp.write('\n\n')
fp.write('Log         = log.$(Process)\n')
fp.write('Output      = out.$(Process)\n')
fp.write('Error       = err.$(Process)\n\n')
fp.write('Queue '+str(number_of_jobs))
fp.close()

print 'Script for BankEfficiency code associated to condor'
print '1 ................................................. BankEfficiency_condor_submit_job created'
print ' a - use the following command to use it : \'condor_submit BankEfficiency_condor_submit_file\''
print ' b - Once the condor script is finished, you\'ll get '+str(number_of_jobs)+' files. '
print '     Merge them with a  command like ''cat out.* > Trigger.dat''\n'
command = path + executable_name + space + arguments +' --print-prototype '
print '\nNow we try to create the xml prototype which might be needed. (same command as in the condor'
print 'script but with the option \'--print-prototype\'):\n'
print ' \t'+command +'\n'
os.system(command)

print '2 ................................................................... xml  prototype created'
print '  its name is BE_Proto.xml'
print '3 Once the condor script is finished and the ascii file \'Trigger.dat \' with all the results'
print 'available you can create the xml file if you wish. Use the  BEAscii2Xml file provided in the'
print 'lalapps/src//findchirp/ package'


