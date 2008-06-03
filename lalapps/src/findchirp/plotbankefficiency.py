#!/usr/bin/python

__Id__ = "$Id$"
__author__ = "Thomas Cokelaer, Thomas.Cokelaer@astro.cf.ac.uk "
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]
__name__ = "plotbankefficiency"
__title__ = "Figure of merits for BankEfficiency results."



  
import getopt, sys, os, re, glob, exceptions, dircache, string
from types    import *
from optparse import *
from matplotlib.ticker import FormatStrFormatter
from numpy import *
from pylab import *
from matplotlib import *
from inspiraltools  import *
#import scipy



#######################################################################
def mysavefig(opts, title):  
  newtitle = title.replace('plotbankefficiency', \
                        'plotbankefficiency_'+opts.user_tag)
  if opts.verbose:
    print '---saving picture ' + newtitle
  
  savefig(newtitle)
##############################################################################
# parse options and arguments

def ParseParameters():
  usage = """Usage: %prog [options] [trigs1 missed1 trigs2 missed2]

  Generate found and missed trig plots

  Example: plotinspmissed --time-dist --show-plot found.xml missed.xml
  """
  parser = OptionParser( usage=usage, \
      version="%prog CVS $Id$ \n" \
      + "$Name$\n" )
  parser.add_option("-s","--show-plot",action="store_true",default=False,\
      help="display the figures on the terminal", dest="show_plot" )
  parser.add_option("-v","--verbose",action="store_true",\
      default=False, help="print information" )
  parser.add_option("-g","--glob",action="store", type="string",\
      default=None, help="print information" )
  parser.add_option("-u","--user-tag",action="store", type="string",\
      default="", help="set a user tag for the figure's filename" )
  parser.add_option("-a","--ambiguity",action="store", type="string",\
      default="", help="the ambiguity function output file from BankEfficiency")
  
  
  (opts,args) = parser.parse_args()

  
  return opts,args
  
command_line = sys.argv[1:]
(opts,args) = ParseParameters()

# Change to Agg back-end if show() will not be called
if not opts.show_plot:
  import matplotlib
  matplotlib.use('Agg')
from pylab import *

rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
rc('text', usetex=True)


# read the xml file to extract the tables
xml = ReadXMLFiles(opts)
results, bank, params,values = xml.getTables()
# now, we can easily extract any parameters from params_table
signal = xml.getvalue("--signal ")

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------  
#the figures -------------------------------------------------------------------
plotting = Plotting(opts.verbose)
plotting.settitle(signal +' injections')
#------------------------------------------------------ overlap versus parameter
# the parameters
params = ['totalmass_sim','eta_sim','chirpmass_sim','inclination_sim',\
          'polarisation_sim','phase_sim']
labels = ['Total Mass ($M_\odot$)','eta','chirp mass ($M_\odot$)',\
          'inclination','polarisation','phase']
for param,label in zip(params,labels):
  # plot the parmeter verus SNR, set the label and title
  plotting.plot(results[param],results['snr'])
  xlabel(r'Injected '+label, size='x-large')
  ylabel(r'Overlap (\%)', size='x-large')
  mysavefig(opts,'plotbankefficiency_snr_versus_'+param+'.png')



#------------------------------------- position of the injections in t0/t3 plane
plotting.plot(results['tau0_sim'],results['tau3_sim'], \
              marker='bo',linewidth=1)
plotting.hold = True
plotting.plot(results['tau0'],results['tau3'], \
              markersize=10, marker='rx', linewidth=0)
try:
  plot(bank['tau0'],bank['tau3'],  markersize=10, marker='x', linewidth=0)
except:
  pass
plotting.hold = False
xlabel(r'$\tau_0$ (seconds) ', size='x-large')
ylabel(r'$\tau_3$ (seconds)', size='x-large')
legend(['Simulated injection','template''s position'])
title('Position of the injections')
mysavefig(opts,'plotbankefficiency_bank.png')

#----------------------------------------------------------SNR histogram and fit
for this,label in zip(\
    ['snr','mass1_sim', 'mass2_sim', 'ecc_sim', 'totalmass_sim',
     'polarisation_sim','inclination_sim'],\
    ['SNR', 'Mass1 (simulation)', 'Mass2(simulation)', \
     'Eccentricity(simulation)','total mass','polarisaton','inclination']):    
  plotting.plot_histogram_and_fit(results[this],100,fit=False)
  xlabel(label, size='x-large')
  ylabel(r'\#', size='x-large')
  mysavefig(opts,'plotbankefficiency_hist_'+this+'.png')

# --------------------------------- scatter plot of SNR versus eccentricty and M
if signal=='Eccentricity':
  # eccentricity, SNr and total mass
  plotting.scatter(results['ecc_sim'],results['totalmass_sim'], \
                   results['snr'], markersize=30, alpha=1,vmin=0,vmax=1)
  xlabel(r'Eccentricity')
  ylabel(r'TotalMass ($M_\odot$)')
  mysavefig(opts,'plotbankefficiency_scattersnr_totalMass_versus_ecc.png')

  # the eccentricity, and SNR
  plotting.plot(results['ecc_sim'],results['snr'],'ob')
  xlabel(r'Eccentricity')
  ylabel(r'Overlap(\%)')
  #pylab.axis([0,0.4,0.7,1])
  mysavefig(opts,'plotbankefficiency_snr_versus_ecc.png')


# --------------------------------------------------------------- accuracies
for param in ['totalmass', 'eta', 'chirpmass','tau0','phase','ecc']:
  plotting.plot(results[param]-results[param+'_sim'], results['snr'],  \
              markersize=10, marker='o', linewidth=0)
  xlabel(r'$\Delta$  ' + param)
  ylabel(r'$\rho$')
  mysavefig(opts,'plotbankefficiency_accuracy_snr_versus_'+param+'.png')


#------------------------------------------------------------------------ PDF
try:
  plotting.vectors2image(results['totalmass_sim'],results['snr'],N=20)
  title(plotting.title+': probability density function of the overlaps')
  mysavefig(opts,'plotbankefficiency_PDF_snr_totalmass_sim.png')
except:
  pass

# ----------------------------------------------------related to the fast option.
plotting.plot(results['totalmass_sim'],results['nfast'])
title(plotting.title+': number of iterations to obtain the overlap computation')
mysavefig(opts,'plotbankefficiency_totalmass_sim_fastoption.png')



# ----------------------------------------------------accuracy of the chirp mass
data = (results['chirpmass_sim']-results['chirpmass'])/results['chirpmass_sim']
plotting.plot_histogram_and_fit(data,50,False)
xlabel(r'$\Delta \mathcal{M}$')
ylabel(r'\#')
mysavefig(opts,'plotbankefficiency_accuracy_chirpmass.png')



# surf as a test function
try:
  x = results['tau0_sim']-results['tau0']
  y = results['tau3_sim']-results['tau3']
  z = results['snr']
  plotting.surf(x,y,z,20,20)
  colorbar()
  mysavefig(opts,'test.png')
except: pass




if opts.show_plot:
  show()





