#!/usr/bin/python 
"""
python script to create a condor_script and xml prototype  
arguments might be changed directly in this file.
"""

__author__ = 'Thomas Cokelaer <Thomas.Cokelaer@astro.cf.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'


import sys
import os
import optparse
import locale
import math
from optparse import OptionParser


path  		= '/home2/spxtc/lscsoft/lalapps/src/findchirp/'
executable_name = 'lalapps_BankEfficiency'


def set_predefined_search_parameter(BE):
	if BE['search']=='BNS':
	    BE['bank-mass-range'] = '1 3'
	    BE['signal-mass-range'] = '1 3'
	    BE['bank-ffinal']= BE['sampling']/2 - 1
    
	elif BE['search']=='BBH':
	    BE['bank-mass-range'] = '3 30'
	    BE['signal-mass-range'] = '3 30'
	    BE['sampling'] = 2048
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	
	elif BE['search']=='BHNS':
	    BE['bank-mass-range'] = '1 33'
	    BE['signal-mass-range'] = '1 30'
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	   
	elif BE['search']=='PBH':
	    BE['bank-mass-range'] = '.3 1'
	    BE['signal-mass-range'] = '.3 1'
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	else:
	    BE['bank-mass-range'] = '1 3'
	    BE['signal-mass-range'] = '1 3'
	    BE['sampling'] = 4096
	    BE['bank-ffinal'] = BE['sampling']/2 - 1
	return BE

def create_condor_file(BE, arguments):
	fp = open('BankEfficiency.sub','w');
	fp.write('Executable   = ' + path + executable_name +'\n')
	fp.write('Universe     = vanilla\n')
	#fp.write('Universe     = vanilla\n')
	#fp.write('Environment  = LD_LIBRARY_PATH=/software/geopptools/lib\n')
	#fp.write('Requirements = Memory >=128 && OpSys == "LINUX" && FileSystemDomain == "explorer" && UidDomain == "explorer"\n')
	#fp.write('+MaxHours =40\n\n')
	fp.write('Arguments = '+ '--seed $(macroseed) ' + arguments)
	fp.write('\n\n')
	fp.write('priority = 10\n')
	
	tag = str(BE['noise-model'])+'_'+str(BE['fl'])+'_'+ str(BE['search']) +'_'+str(BE['signal'])+'_'+str(BE['signal-order'])+'_'+str(BE['template'])+'_'+str(BE['template-order'])+'_'+str(BE['sampling'])+'_'+str(BE['mm'])+'.$(macroseed)\n'
	msg = 'log = ./log/log_'+tag
	fp.write(msg)
	msg = 'output = out_'+tag
	fp.write(msg)
	msg = 'error = ./log/err_'+tag
	fp.write(msg)
	fp.write('notification = never\n')

	fp.write('Queue 1')
	fp.close()

def create_bank(arguments):
	os.system('rm -f BE_Bank.dat BE_Bank.xml')
	print '###'
	print ' We are creating the template bank for sanity check. Please wait'
	fp =open('BankEfficiency_createbank','w');
	fp.write( path + executable_name + arguments+' --n 1 --faithfulness --print-bank 1> out 2>err'+'\n')
	fp.close()
	os.system('chmod 755 BankEfficiency_createbank')
	a=os.system('./BankEfficiency_createbank')
	
	if a==0:
		print '->done (your parameters seems correct). See BE_Bank.xml file.'
	else:
		print '->failed (your parameters seems correct)'
		quit

def create_dag_file(njobs):
	fp=open('BankEfficiency.dag', 'w')
	for id in range(1,njobs+1,1):
		fp.write('JOB '+str(id)+' BankEfficiency.sub'+'\n')
	 	fp.write('VARS '+str(id)+' macroseed="'+str(id)+'"\n')
	fp.close()


def check_executable():
	try:
		print 'check that the excutable is present ...'
		f = open(path+executable_name, 'r')
		f.close()
	except:
		print 'Can not find ' +path + executable_name
		sys.exit()
	print 'lalapps_bankefficiency found. Going ahead'
	print '--'
	



def main():

    parser = OptionParser()
    
    parser.add_option("", "--noise-model", 
		dest='noise_model',default='LIGOI',metavar='NOISEMODEL',
		help=" <VIRGO, GEO, LIGOI, LIGOA>") 
    parser.add_option("", "--search",
		dest='search',default='BNS',
		help=" <BNS, BBH, PBH , BHNS>")
    parser.add_option("","--signal-mass-range",
		default='1 3', dest='signal_mass_range',
		help="min and max individual mass in solar mass." )
    parser.add_option("","--bank-mass-range",
		default='1 3', dest='bank_mass_range',
		help="min and max individual mass in solar mass." )
    parser.add_option("","--signal",
		default='EOB', dest='signal',type='string',
		help="approximant of the injection (EOB, TaylorT1, ...)." )
    parser.add_option("","--template",
		default='EOB', dest='template',type='string',
	 	help="approximant of the injection (EOB, TaylorT1, ...)." )
    parser.add_option("","--template-order",
		default=4, dest='template_order',type='int',
		help="PN order of the template." )
    parser.add_option("","--signal-order",
		default=4,dest='signal_order', type='int',
		help="PN order of the signal." )
    parser.add_option("","--minimal-match",
		default=0.95, dest='minimal_match', 
		help="minimal match." )
    parser.add_option("","--sampling",
		dest='sampling', default=4096, type='float',
		help="sampling frequency" )
    parser.add_option("","--bank-ffinal",
		dest='bank_ffinal', default=2047, type='float',
		help="upper frequency to be used" )
    parser.add_option("-n","--ntrial",
		dest='ntrial', default=10000, type='int',
		help="number of trial." )
    parser.add_option("","--njobs",
		dest='njobs', default=100, type='int',
		help="number of jobs." )
    parser.add_option("","--bank-grid-spacing",
		dest='bank_grid_spacing', default='Hexagonal', 
	 	help="type of template bank placement : Hexagonal, SquareNotOriented, HexagonalNotOriented" )
    parser.add_option("","--fl",
		dest='fl',  type='int',default=-1, 
		help="lower cut off frequency" )
    parser.add_option("","--max-total-mass",
		dest='max_total_mass', default=-1, type='float',
		help="max total mass (injection)" )
    parser.add_option("","--fast-simulation",
		default="false",
		dest='fast_simulation', 
		help="fast simulation" )


    (options, args) = parser.parse_args()
    BE={} 
    BE['search'] = options.search 
    BE['noise-model']=options.noise_model
    BE['fl']=options.fl
    BE['bank-mass-range']=options.bank_mass_range
    BE['signal-mass-range']=options.signal_mass_range
    BE['signal']=options.signal
    BE['template']=options.template
    BE['signal-order']=options.signal_order
    BE['template-order']=options.template_order
    BE['mm']=options.minimal_match
    BE['sampling']=options.sampling
    BE['bank-grid-spacing']=options.bank_grid_spacing
    #derived options set to fixed value for the time being.
    if options.bank_ffinal > options.sampling/2-1:
	BE['bank-ffinal'] = options.sampling/2-1
    else:
	BE['bank-ffinal'] = options.bank_ffinal
    
    #[some other default values]
    others = ' --print-result-xml  --debug 33 --print-bank '
    if options.fast_simulation==True:
  	others = others + ' --fast-simulation '
    arguments  = others

    #depending on the "search" value, we set some extra default values
    BE = set_predefined_search_parameter(BE)
    #check that the executable is present
    check_executable()

    if BE['fl']==-1:
        if BE['noise-model']=='VIRGO':
            BE['fl'] = 20        
        elif BE['noise-model']=='GEO':
            BE['fl']=40
        elif BE['noise-model']=='LIGOI':
            BE['fl']=40
        elif BE['noise-model']=='LIGOA':
            BE['fl']=20

    # compute the number of trial per node   
    nCondor = math.ceil(options.ntrial/options.njobs)
    BE['ntrial'] = nCondor

    for arg in  BE:
	if arg!='search':
		arguments = arguments +  ' --'+arg+' '+ str(BE[arg])

    if options.max_total_mass > 0:
	arguments = arguments + ' --max-total-mass ' + str(option.max_total_mass)

    # print some information on the screen
    print """
	The condor script will use the following arguments 
	-------------------------------------------
    """

    for arg in BE:
	print '    ' + arg + ' = ' +str(BE[arg])
    print 'The number of simulation requested is '+str(BE['ntrial'])
    print 'They will be split into '+ str(options.njobs)+' jobs'

    # create the condor file using the input parameter stored in BE
    create_condor_file(BE, arguments)
    create_dag_file(options.njobs)

    # we create only the bank. This is mainly to test if the 
    # input parameters are correct
    create_bank(arguments)

    print 'Generating the prototype xml file for merging condor job'
    command = path + executable_name + ' ' + arguments +' --print-prototype 1>out 2>err'
    os.system(command)
    print '->done'
    
    print """ In order to start the job, type
    	-->	condor_submit BankEfficiency.dag
    	which reads BankEfficiency.sub. Once the dag is finished, 
    	you\'ll get '+str(options.njobs)+' files. 
    	concatenate the job using 
    	cat out.* > Trigger.dat
    	cp BE_Proto.xml and TMPLTBANK.xml into your directory and type
    	./lalapps_bankefficiency --ascii2xml """




if __name__ == "__main__":
    main()
