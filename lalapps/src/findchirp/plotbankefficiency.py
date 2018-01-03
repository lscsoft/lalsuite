#!/usr/bin/python
__Id__ = "$Id$"
__author__ = "Thomas Cokelaer, Thomas.Cokelaer@astro.cf.ac.uk "
__version__ = "$Revision$"
__date__ = "$Date$"
__name__ = "plotbankefficiency"
__title__ = "Figure of merits for BankEfficiency results."
  
import getopt, sys, os, re, glob, exceptions, dircache, string
from optparse import *
from matplotlib.ticker import FormatStrFormatter
from numpy import *
from pylab import *
from matplotlib import rc
from lalapps.inspiraltools import *

#this is for the colorbar labelsize. Do not know why it works but it works!
rc('ytick',labelsize=15)


def setTickSize(fontsize):
  """
  An alias fonction to set the fontsize of the ticks 
  """
  ax = gca()
  for tick in ax.xaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
  for tick in ax.yaxis.get_major_ticks():
    tick.label1.set_fontsize(fontsize)
    

def fastoptionRelated(results):
  if opts.verbose:
    print 'Entering fastoptionRelated ---------------------------------------- '
  # ----------------------------------------------------related to the fast option.
  name = __name__+'_totalmass_sim_fastoption'
  try:
    plotting.plot(results['totalmass_sim'],results['nfast'])
    plotting.hold = True
    plotting.plot(results['totalmass_sim'],results['nfast_max'], marker='or')
    plotting.hold = False
    title(plotting.title+': number of iterations performed to obtain the overlap')
    mysavefig(opts,name + '.png')
  
    plotting.plot(results['totalmass_sim'],results['nfast']/results['nfast_max'])
    mysavefig(opts,name + '_ratio.png')
  except:pass  
  try:
    name = __name__+'_scatter_ecc_sim_totalmass_sim_fastoption_ratio'
    plotting.scatter(results['ecc_sim'],results['totalmass_sim'],\
      results['nfast']/(results['nfast_max']+1), markersize=30, alpha=1)
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'TotalMass ($M_\odot$)',fontsize=opts.fontsize)
    title('Ratio between the total number of filter and the number of the filter that gives the best SNR')
    mysavefig(opts, name+'.png')
  except:pass
  
  
  try:
    name = __name__+'_scatter_ecc_sim_totalmass_sim'
    plotting.scatter(results['ecc_sim'],results['totalmass_sim'],\
      results['nfast_max'], markersize=30, alpha=1)
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'TotalMass ($M_\odot$)',fontsize=opts.fontsize)
    mysavefig(opts, name+'_nfast_max.png')
  except:pass
  try:
    plotting.scatter(results['ecc_sim'],results['totalmass_sim'],\
      results['nfast'], markersize=30, alpha=1)
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'TotalMass ($M_\odot$)',fontsize=opts.fontsize)
    mysavefig(opts, name+'_nfast.png')
  except:pass
  
    
  

def EMatchRelated(results):
  if opts.verbose:
    print 'Entering EMatchRelated -------------------------------------------- '
  
  name = __name__+'_snr_versus_EMatch'
  try:
    clf()
    plotting.plot(results['ematch'],results['snr'])
    xlabel('EMatch',fontsize=opts.fontsize)
    ylabel(r'Overlap (\%)', fontsize=opts.fontsize)
    plotting.hold = True
    plotting.plot([min(results['ematch']), max(results['ematch'])],[ematch,ematch], marker='r',linewidth=2)
    mysavefig(opts,name + '.png')
    plotting.hold = False
  except:
    print 'Problem inside plotting.plot '+ name
    pass


def accuraciesRelated(results):
  if opts.verbose:
    print 'Entering accuraciesRelated ---------------------------------------- '
  
  # --------------------------------------------------- accuracies
  for param in ['totalmass', 'eta', 'chirpmass','tau0','tau3','phase','ecc','ffinal']:
    name = __name__+'_accuracy_snr_versus_'+param
    try:
      plotting.plot(results[param]-results[param+'_sim'], results['snr'],  \
              markersize=10, marker='o', linewidth=0)
      xlabel(r'$\Delta$  ' + param,fontsize=opts.fontsize)
      ylabel(r'$\rho$',fontsize=opts.fontsize)
      mysavefig(opts, name+'.png')
    except:
      print 'Problem inside plot '+name
      pass


def HistogramRelated(results):
  if opts.verbose:
    print 'Entering HistogramRelated ----------------------------------------- '
    
  params = ['snr','mass1_sim', 'mass2_sim', 'ecc_sim', 'totalmass_sim',
     'polarisation_sim','inclination_sim']
  labels = ['SNR', 'Mass1 (simulation)', 'Mass2(simulation)', \
     'Eccentricity(simulation)','total mass','polarisaton','inclination']
  for this,label in zip(params, labels):
    name = __name__+'_hist_'+this
    try:
      plotting.plot_histogram_and_fit(results[this],100,fit=False)
      xlabel(label, fontsize=opts.fontsize)
      ylabel(r'\#', fontsize=opts.fontsize)
      mysavefig(opts,name+'.png')
    except:
      print 'Problem inside plot_histogram_and_fit' + name
      pass




def BankRelated(results, bank):
  if opts.verbose:
    print 'Entering BankRelated ---------------------------------------------- '
  
  #------------------------------------- position of the injections in t0/t3 plane
  plotting.plot(results['tau0_sim'],results['tau3_sim'], \
              marker='bo',linewidth=1)
  plotting.hold = True
  plotting.plot(results['tau0'],results['tau3'], \
              markersize=10, marker='rx', linewidth=0)
  try:
    plot(bank['tau0'],bank['tau3'],  markersize=15, marker='x', linewidth=0)
  except:
    pass
  plotting.hold = False
  xlabel(r'$\tau_0$ (seconds) ', fontsize=opts.fontsize)
  ylabel(r'$\tau_3$ (seconds)', fontsize=opts.fontsize)
  legend(['Simulated injection','template''s position'])
  title('Position of the injections')
  mysavefig(opts,'plotbankefficiency_bank.png')

def SNRRelated(results):
  if opts.verbose:
    print 'Entering SNRRelated ----------------------------------------------- '
  
  params = ['totalmass_sim','eta_sim','chirpmass_sim','inclination_sim',\
          'polarisation_sim','phase_sim']
  labels = ['Total Mass ($M_\odot$)','eta','chirp mass ($M_\odot$)',\
          'inclination','polarisation','phase']
  for param,label in zip(params,labels):
    # plot the parameter versus SNR, set the label and title
    name = __name__+'_snr_versus_'+param
    try:
      plotting.plot(results[param],results['snr'])
      xlabel(r'Injected '+label,fontsize=opts.fontsize)
      ylabel(r'Overlap (\%)',fontsize=opts.fontsize)
      mysavefig(opts, name + '.png')
    except:
      print 'Problem inside plotting.plot '+ name
      pass



def comparisonRelated(results, results2):
  if opts.verbose:
    print 'Entering comparisonRelated ---------------------------------------- '
  try:
    clf()
    plotting.plot(results['totalmass_sim'], results['snr']/results2['snr'])
    xlabel('total mass',fontsize=opts.fontsize)
    ylabel('SNR in file 1 /  SNR in file 2',fontsize=opts.fontsize)
    mysavefig(opts, 'plotbankefficiency_compare_totalmass_versus_snr_ratio.png')
  except: 
    print 'Problem in comparisonRelated. Maybe the two files do not have the same size.'
    pass
  try:
    clf()
    plotting.plot(results['totalmass_sim'], results['snr'])
    plotting.hold = True
    plotting.plot(results2['totalmass_sim'], results2['snr'], marker='rx')
    plotting.hold = False
    xlabel('total mass',fontsize=opts.fontsize)
    ylabel('SNR in file 1 /  SNR in file 2',fontsize=opts.fontsize)
    mysavefig(opts, 'plotbankefficiency_compare_totalmass_versus_snr.png')
  except: 
    print 'Problem in comparisonRelated. Maybe the two files do not have the same size.'
    pass


def bcvRelated(results):
  if opts.verbose:
    print 'BCV related plotting ---------------------------------------------- '
  try:
    plotting.plot(results['totalmass_sim'],results['alpha_f'])
    xlabel('total mass',fontsize=opts.fontsize)
    ylabel('alpha\_F',fontsize=opts.fontsize)
    mysavefig(opts, 'plotbankefficiency_alphaf_vs_totalmass.png')
  except:pass



def eccentricityRelated(results, results2, nbin=40):
  if opts.verbose:
    print 'Eccentricity related plotting ------------------------------------ '

  totMass = results['totalmass_sim']
  snr = results['snr'] 
  ecc   = xml.getvalue('--signal-eccentricity-range')
  eccRange = [float(x) for x in ecc.split(" ")]
  eccMin = eccRange[0]
  eccMax = eccRange[1]*1.01
 
  # Scatter plots of the eccentricity versus SNR and total mass
  name = __name__ + '_scatter_snr_versus_totalmass_ecc'
  try:
    plotting.scatter(results['ecc_sim'],totMass, snr,  markersize=30, alpha=1,\
                     vmin=0,vmax=1)
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'TotalMass ($M_\odot$)',fontsize=opts.fontsize)
    mysavefig(opts,name+'.png')
  except:
    print 'Problem inside scatter1 ' +name
    pass

  if opts.compare_with is not None:
    snr2 = results2['snr']
    totMass2 = results2['totalmass_sim']
    name = __name__ + '_compare_scatter_snr_versus_totalmass_ecc'    
    try:
      plotting.scatter(results['ecc_sim'],totMass, snr/snr2,  \
          markersize=30, alpha=1,vmin=0.8,vmax=1.2)
      xlabel(r'Eccentricity',fontsize=opts.fontsize)
      ylabel(r'TotalMass ($M_\odot$)',fontsize=opts.fontsize)
      mysavefig(opts, name+'.png')
    except:
      print 'Problem inside scatter2 ' +name,
      pass
  if opts.compare_with is not None:
    name = __name__ + '_compare_snr_versus_totalmass_ecc'    
    try:
      plotting.plot(results['ecc_sim'],snr)
      plotting.hold = True
      plotting.plot(results2['ecc_sim'],snr2, marker='rx')
      plotting.hold = False
      xlabel(r'Eccentricity')
      ylabel(r'Overlap (\%)')
      mysavefig(opts,name+'.png')
    except:
      print 'Problem inside plot ' +name
      pass


  # Simple 2D plot of the SNR versus eccentricity 
  name = __name__+'_snr_versus_ecc'
  try:
    plotting.plot(results['ecc_sim'],snr,'ob')
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'Overlap(\%)',fontsize=opts.fontsize)
    #pylab.axis([0,0.4,0.7,1])
    mysavefig(opts,name+'.png')
  except:
    print 'Problem inside plot ' +name
    pass

  # contour plots of SNR versus eccentricity (at 2*fl/3) and SNR
  if float(xml.getvalue('--m1'))!=-1:
    return

  name = __name__ + '_surf_snr_versus_totalmass_ecc'

  try:
    plotting.surf(results['ecc_sim'], totMass, snr,\
                  xbin=nbin,ybin=nbin,vmin=0.5,vmax=1)
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'Total mass $(M_\odot)$',fontsize=opts.fontsize)
    mysavefig(opts,name+'.png')
  except:
    print 'Problem inside surf ' + name
    pass
  
  # create 4 plots for chirpmass and totalMass
  # and for ecc_sim and ecc_sim_fl (at 2*fl/3)
  for ecc in ['ecc_sim','ecc_sim_fl']:
    for this_param,label in zip(['totalmass_sim','chirpmass_sim'],\
        [r'Total mass $(M_\odot)$',r'$\mathcal{M} (M_\odot)$']):
#      try:
        name = __name__ + '_contour_snr_versus_'+this_param+'_'+ecc
        h = plotting.contourf(results[ecc], results[this_param], snr, \
            xbin=nbin,ybin=nbin, xmin=eccMin, xmax=eccMax, vmin=0.5,vmax=1,\
            colorbar_title=r'Overlap(\%)',fontsize=opts.fontsize)
        xlabel(r'Eccentricity',fontsize=opts.fontsize)
        ylabel(label,fontsize=opts.fontsize)
        setTickSize(opts.fontsize/1.2)
        
        mysavefig(opts,name +'.png')
        
        name = __name__ + '_contour_detectability_versus_'+this_param+'_'+ecc
        h = plotting.contourf(results[ecc], results[this_param], (snr)**(3.),\
                      xbin=nbin, ybin=nbin, xmin=eccMin, xmax=eccMax, vmin=0.5,vmax=1,colormap='cube',\
                      colorbar_title=r'Overlap(\%)',fontsize=opts.fontsize)
        xlabel(r'Eccentricity',fontsize=opts.fontsize)
        ylabel(label,fontsize=opts.fontsize)
        setTickSize(opts.fontsize/1.2)
        
        mysavefig(opts,name +'.png')        
 #     except: 
 #       pass
  
    
  
  try:
    name = __name__ +'_eccentricity_ratio_tau0'
    clf()
    plotting.scatter(results['ecc_sim'],results['totalmass_sim'],\
                     results['tau0_sim']/results['tau0'],markersize=30, alpha=1)
    xlabel(r'Eccentricity')
    ylabel(r'Total mass $(M_\odot)$')
    title(r'Simulated $\tau_0$ divided by estimated $\tau_0$')
    mysavefig(opts,name+'.png')  
  except: 
    print 'Problem in '+ name
    pass

  try:
    name = __name__ +'_eccentricity_ratio_tau3'
    clf()
    plotting.scatter(results['ecc_sim'],results['totalmass_sim'],\
                     results['tau3_sim']/results['tau3'],markersize=30, alpha=1)
    xlabel(r'Eccentricity',fontsize=opts.fontsize)
    ylabel(r'Total mass $(M_\odot)$',fontsize=opts.fontsize)
    title(r'Simulated $\tau_3$ divided by estimated $\tau_3$')
    mysavefig(opts,name+'.png')  
  except:
    print 'Problem in '+ name
    pass

  
  

#######################################################################
def mysavefig(opts, title):  
  newtitle = title.replace('plotbankefficiency', \
                        'plotbankefficiency_'+opts.user_tag)
  if opts.verbose:
    print 'saving picture ' + newtitle
  
  savefig(newtitle)
##############################################################################
# parse options and arguments

def ParseParameters():
  usage = """Usage: %prog [options]

  Generate figures of merits using the XML output of BankEfficieny

  Example: python plotbankefficiency.py --glob 'BankEfficiency-Results.xml' --verbose --user-tag tag1
  """
  parser = OptionParser( usage=usage, \
      version="%prog CVS $Id$ \n" \
      + "$Name$\n" )
  parser.add_option("-s","--show-plot",action="store_true",default=False,\
      help="display the figures on the terminal", dest="show_plot" )
  parser.add_option("-v","--verbose",action="store_true",\
      default=False, help="print information" )
  parser.add_option("-g","--glob",action="store", type="string",\
      default=None, help="name of the file to investigate" )
  parser.add_option("-c","--compare-with",action="store", type="string",\
      default=None, help="a file to compare " )
  parser.add_option("-u","--user-tag",action="store", type="string",\
      default="", help="set a user tag for the figure's filename" )
  parser.add_option("-a","--ambiguity",action="store", type="string",\
      default="", help="the ambiguity function output file from BankEfficiency")
  parser.add_option("-n","--n",action="store", type="float",\
      default="1e16", help="number of simulation to look at")
  parser.add_option("--skip-accuracy",action="store_true",\
      default=False, help="skip plots related to parameter accuracy")
  parser.add_option("--skip-fast",action="store_true", \
      default=False, help="skip plots related to fast option")
  parser.add_option("--skip-eccentricity",action="store_true",\
      default=False, help="skip plots related to eccentricity")
  parser.add_option("--skip-ematch",action="store_true",\
      default=False, help="skip plots related to ematch")
  parser.add_option("--skip-bcv",action="store_true",\
      default=False, help="skip plots related to BCV")
  parser.add_option("--skip-snr",action="store_true",\
      default=False, help="skip plots related to 2D plots with SNR")
  parser.add_option("--skip-bank",action="store_true",\
      default=False, help="skip plots related to the bank")
  parser.add_option("--skip-histogram",action="store_true",\
      default=False, help="skip plots related to histograms")
  parser.add_option("--nbin",action="store",\
      default=40, help="set binning of the contour plots")
  parser.add_option("--fontsize",action="store",type="float",\
      default=20, help="set fontsize of the labels")
  
  
  
  (opts,args) = parser.parse_args()

  
  return opts,args
  
  

if __name__ == 'plotbankefficiency': 
  
  command_line = sys.argv[1:]
  (opts,args) = ParseParameters()

  # Change to Agg back-end if show() will not be called
  if not opts.show_plot:
    import matplotlib
    matplotlib.use('Agg')
  from pylab import *

  rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
  rc('text', usetex=True)
  
  # -------------------------------------------------------- THE DATA
  # Read the XML file
  xml = ReadXMLFiles(opts)
  # Read the tables (bankefficiency, sngl-inspiral, process_params)
  results, bank, params,values = xml.getTables()
  # extract some parameters
  signal   = xml.getvalue("--signal ")
  template = xml.getvalue("--template ")
  ematch   = xml.getvalue('--e-match')
  ecc   = xml.getvalue('--signal-eccentricity-range')
  eccRange = [float(x) for x in ecc.split(" ")]
# ---------------------------------------------------------- THE FILE TO COMPARE
  if opts.compare_with is not None:
    opts.glob = opts.compare_with
    xml2      = ReadXMLFiles(opts)
    results2,bank2,param2,values2 = xml2.getTables()
  else: 
    results2 = None

  # ----------------------------------------------------------- PLOTTING

  plotting = Plotting(opts.verbose)
  plotting.settitle(signal +' injections')
  #------------------------------------------------------ overlap versus parameter

  if opts.skip_snr is False:
    try: SNRRelated(results)
    except: print 'Problem in SNRRelated)'

  if opts.skip_bank is False:
    try: BankRelated(results,bank)
    except: print 'Problem in BankRelated'
  
  if opts.skip_histogram is False:  
    try: HistogramRelated(results)
    except: print 'Problem in HistogramRelated'

  if opts.skip_eccentricity is False:
    if signal=='Eccentricity':
    #try:
      eccentricityRelated(results,results2,nbin=opts.nbin)
    #except:print 'Problem in eccentricityRelated'

  if opts.skip_bcv is False:
    if template=='BCV':
      bcvRelated(results)

# --------------------------------------------------- compare-with a second file
  if results2 is not None:
    try: comparisonRelated(results, results2)
    except: print 'Problem in comparison related'
  
  if opts.skip_accuracy is False:
    try: accuraciesRelated(results)
    except: print 'Problem in accuraciesRelated'

  if opts.skip_fast is False:
    try: fastoptionRelated(results)
    except: print 'Problem in fast option related'

  if opts.skip_ematch is False:
    try: EMatchRelated(results)
    except:print 'Problem in EMatchRelated'


  if opts.show_plot:
    show()





