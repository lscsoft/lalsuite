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
path  		= './'
executable_name = 'lalapps_BankEfficiency'
# 
# injection parameters for TD
min_mass       	=  3  	   	# miminal individaul mass to inject
max_mass       	=  80  	# maximal individua lmass to inject
maxTotalMass    = -1
# bank parameters for TD 
min_massb      	=  3  		# maximal individua lmass to inject
max_massb      	=  80  	# maximal individua lmass to inject
#injection parameters for BCV
signal_psi0_min	= 1000
signal_psi0_max	= 550000
signal_psi3_min	= -5000
signal_psi3_max	= -1000
#bank parametesr for BCV 
bank_psi0_min	= 10000
bank_psi0_max	= 1500000
bank_psi3_min	= -7000
bank_psi3_max	= -10
alpha_bank	= 0.01
fend_min	=  -2  
fend_max	=  6 
number_fcut	=  3 
# common bank parameters 
minimal_match 	 = 0.95     	# minimal match in bank creation
lower_frequency  = 40		# lower freqency cutoff. identical for both signal and templates for the time being
freq_moment_bank = 2047
sampling = 4096
#injection
signal 		= 'TaylorT3'
signal_order 	=  4
# template bank
template 	= 'TaylorT3'
template_order  =  4

noisemodel	= 'LIGOI'
simulation_type = 'SignalOnly'       	# 0 for signal only 1 for noise only and 2 for noise+signal
signal_amplitude = 8 

#condor parameters 
number_of_jobs 	= 500
number_of_trials_per_job = 200

#real data parameters 
channel = "H1:LSC-AS_Q"
startTime = 752058105
#startTime = 755512721
#startTime = 757575852


others = '--print-result-xml  --bank-grid-spacing Hexagonal  --debug 33 --print-bank --bank-inside-polygon 0 --data-checkpoint --fast-simulation'


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




arguments = ' --n	' 	  + str(number_of_trials_per_job) \
         +' --max-total-mass '  +str(maxTotalMass) \
         +' --signal-mass-range ' + str(min_mass)+ space +str(max_mass) \
         +' --bank-mass-range '	  + str(min_massb)+ space +str(max_massb) \
       	 +' --mm ' 		  + str(minimal_match) \
         +' --fl ' 		  + str(lower_frequency) \
         +' --signal ' 		  + str(signal) \
         +' --signal-order ' 	  + str(signal_order) \
         +' --template ' 	  + str(template) \
         +' --template-order ' 	  + str(template_order) \
         +' --simulation-type '	  + str(simulation_type) \
         +' --signal-amplitude '  + str(signal_amplitude) \
         +' --bank-alpha ' 	  + str(alpha_bank) \
	 +' --noise-model ' 	  + noisemodel \
         +' --bank-fcut-range '	  + str(fend_min) + space + str(fend_max) \
         +' --bank-number-fcut '  + str(number_fcut) \
         +' --bank-psi0-range '	  + str(bank_psi0_min) + space + str(bank_psi0_max) \
         +' --signal-psi0-range ' + str(signal_psi0_min) + space + str(signal_psi0_max) \
         +' --bank-psi3-range '	  + str(bank_psi3_min) + space + str(bank_psi3_max)\
         +' --signal-psi3-range ' + str(signal_psi3_min) + space + str(signal_psi3_max)\
         +' --bank-ffinal '       + str(freq_moment_bank)\
         +' --channel '           + str(channel) \
         +' --gps-start-time '    + str(startTime) \
         +' --sampling '          + str(sampling) \
	 +space + others


# create the condor file
fp =open('BankEfficiency_condor_submit_job','w');
fp.write('Executable   = ' + path + executable_name +'\n')
#fp.write('Universe     = vanilla\n')
fp.write('Universe     = standard\n')
#fp.write('Environment  = LD_LIBRARY_PATH=/software/geopptools/lib\n')
#fp.write('Requirements = Memory >=128 && OpSys == "LINUX" && FileSystemDomain == "explorer" && UidDomain == "explorer"\n')
fp.write('+MaxHours =40\n\n')
fp.write('Arguments = '+ '--seed $(Process) ' + arguments)
fp.write('\n\n')
fp.write('priority = 10\n')
fp.write('Log         = log.$(Process)\n')
fp.write('Output      = out.$(Process)\n')
fp.write('Error       = err.$(Process)\n\n')
fp.write('Notification = never\n')
fp.write('Queue '+str(number_of_jobs))
fp.close()


#possibly create a dag file if n is large 
fp =open('BankEfficiency.sub','w');
fp.write('Executable   = ' + path + executable_name +'\n')
#fp.write('Universe     = vanilla\n')
fp.write('Universe     = standard\n')
fp.write('Arguments = '+ '--seed $(macroseed) ' + arguments)
fp.write('\n\n')
fp.write('priority = 10\n')
fp.write('log         = ./log/tmp\n')
fp.write('output      = ./log/out.$(macroseed)\n')
fp.write('error       = ./log/err.$(macroseed)\n\n')
fp.write('notification = never\n')
fp.write('queue 1')
fp.close()
#possibly create a dag file if n is large 
fp=open('BankEfficiency.dag', 'w')
for id in range(1,number_of_jobs+1,1):
 fp.write('JOB '+str(id)+' BankEfficiency.sub'+'\n')
 fp.write('VARS '+str(id)+' macroseed="'+str(id)+'"\n')
fp.close()

# create bank
print '###'
print ' We are creating the template bank for sanity check (could be long with real data)'
fp =open('BankEfficiency_createbank','w');
fp.write( path + executable_name + arguments+' --n 1 --faithfulness --print-bank 1> out 2>err'+'\n')
fp.close()
os.system('chmod 755 BankEfficiency_createbank')
a=os.system('./BankEfficiency_createbank')

if a==0:
	print '... done (your parameters seems correct). See BE_Bank.xml file.'
else:
	print '... failed (your parameters seems correct)'
	quit

print 'Generating the prototype xml file for merging condor job'
command = path + executable_name + space + arguments +' --print-prototype 1>out 2>err'
os.system(command)
print 'done'

print ' creating the condor script ... BankEfficiency_condor_submit_job created'
print '------'
print ' type ''condor_submit BankEfficiency_condor_submit_job'' to launch the job '
print '------'
print ' Once the condor scripts are finished, you\'ll get '+str(number_of_jobs)+' files. '
print ' Merge them with ''cat out.* > Trigger.dat'' and type ''BEAscii2Xml'' which search for BE_Proto.xml and Trigger.dat files'
print ''
