#!/usr/bin/python
__author__ = "Rahul Biswas <rahul@gravity.phys.uwm.edu>, Kipp Cannon <kipp.cannon@ligo.org>, Ruslan Vaulin <vaulin@gravity.phys.uwm.edu>"



import sys
import exceptions
import glob
import exceptions
import matplotlib.numerix.ma as ma
import matplotlib.colors as colors
import glob
import ConfigParser
import string
import fileinput
import linecache
import math

from glue.ligolw import table
from glue.ligolw import lsctables
from glue.ligolw import utils
from glue import segments

from pylal import rate
from pylal import git_version
from pylal import CoincInspiralUtils
from pylal import SnglInspiralUtils
from pylal import ligolw_burca_tailor
from pylal.tools import XLALCalculateEThincaParameter
from optparse import *
from pylab import *

###############################################################################
usage = """
usage: %prog [options]

Script to generate histogram of snrs in various coincidences, for background (time-slides) and injections. 

"""
#
#
#
###############################################################################

# Input options
###############################################################################


parser = OptionParser(
	version = "Name: %%prog\n%s" % git_version.verbose_msg,
	usage
)

parser.add_option("-G","--found-glob",action="store",type="string",\
    default=None, metavar=" FOUND_GLOB", \
    help="GLOB of found trigger/injection files to read" )

parser.add_option("-M","--missed-glob",action="store",type="string",\
    default=None, metavar=" MISS_GLOB", \
    help="GLOB of files containing missed injections to read" )

parser.add_option("-S","--source-file",action="store",type="string",\
    default=None,metavar=" SOURCES",help="full path to source file")

parser.add_option("-I","--injection-glob",action="store",type="string",\
    default=None,metavar=" INJ_GLOB",\
    help="GLOB of files containing astrophysically distributed injections")

parser.add_option("-K","--slide-glob",action="store",type="string",\
    default=None,metavar=" SLIDE_GLOB",\
    help="GlOB of files containing slide triggers" )

parser.add_option("-Z","--zero-glob",action="store",type="string",\
    default=None,metavar=" ZERO_GLOB",\
    help="GlOB of files containing zero lag triggers" )


parser.add_option("-A","--statistic",action="store",default='snr',\
    type="string",\
    help="choice of statistic used in making plots, valid arguments are:"  "snr (DEFAULT), snr_over_chi, s3_snr_chi_stat, effective_snr, "
    "bitten_l, bitten_lsq")

parser.add_option("-P","--param",action="store",default=None,\
    type="string",\
    help="First ifo in double coincidence")


parser.add_option("-R","--ifo",action="store",default=None,\
    type="string",\
    help="Single ifo for which parameters to be plotted")


parser.add_option("-C","--coincs",action="store",default=None,\
    type="string",\
    help="Coincs for doubles")
""

(opts,args) = parser.parse_args()

########################################



###############################################################################
# Checklist to see if all the files have been read in order to proceed further
###############################################################################


# Zero lag  Files
if opts.zero_glob:
  zerolagfiles = []
  zerolagfiles = glob.glob(opts.zero_glob)
  if len(zerolagfiles) < 1:
    print >>sys.stderr, "The glob for " + opts.zero_glob + " returned no files"
    sys.exit(1)

# Time slides Files
if opts.slide_glob:
  slidesfiles = []
  slidesfiles = glob.glob(opts.slide_glob)
  if len(slidesfiles) < 1:
    print >>sys.stderr, "The glob for " + opts.slide_glob + " returned no files"
    sys.exit(1)


# Injection Files
if opts.injection_glob:
  injfiles = []
  injfiles = glob.glob(opts.injection_glob)
  if len(injfiles) < 1:
    print >>sys.stderr, "The glob for " + opts.injection_glob + "returned no files"
    sys.exit(1)


# Found Injections 
if opts.found_glob:
  foundfiles = []
  foundfiles = glob.glob(opts.found_glob)
  if len(foundfiles) < 1:
    print >>sys.stderr, "The glob for " + opts.found_glob + "returned no files "
    sys.exit(1)


# Injection Files
if opts.missed_glob:
  missedfiles = []
  missedfiles = glob.glob(opts.missed_glob)
  if len(missedfiles) < 1:
    print >>sys.stderr, "The glob for " + opts.missed_glob + "returned no files "
    sys.exit(1)


# check that statistic is OK:
if (opts.statistic != 'snr') and (opts.statistic != 'snr_over_chi') \
    and (opts.statistic != 's3_snr_chi_stat') \
    and (opts.statistic != 'effective_snr')\
    and (opts.statistic != 'bitten_lsq')\
    and (opts.statistic != 'bitten_l'):
  print >>sys.stderr, "--statistic must be one of"
  print >>sys.stderr, "(snr|snr_over_chi|s3_snr_chi_stat|effective_snr|bitten_l)"
  sys.exit(1)

###############################################################################
# List of ifos 

ifos = ["H1","H2","L1"]

if opts.coincs == "H1H2":
  ifolist  = ("H1","H2")
elif opts.coincs == "H1L1":
  ifolist = ("H1","L1")
elif opts.coincs == "H2L1":
  ifolist = ("H2","L1")

###############################################################################
# Extract column of interest from the SnglInspiralTable
###############################################################################


###############################################################################


def single_params_func(events, timeslide = 0):

        H1_eff_snr = events.get_effective_snr()
        H2_eff_snr = events.get_effective_snr()
        L1_eff_snr = events.get_effective_snr()
     
        return {
                "H1_eff_snr": H1_eff_snr,
                "H2_eff_snr": H2_eff_snr,
                "L1_eff_snr": L1_eff_snr, 
        }

def double_params_func(events, timeslides = 0):
       
        if hasattr(events, "H1") and hasattr(events, "H2"): 
          H1H2_eff_snr = math.sqrt((events.H1.get_effective_snr())**2 + (events.H2.get_effective_snr())**2)
          return {
                "H1H2_eff_snr": H1H2_eff_snr
          }

         
        elif hasattr(events, "H1") and hasattr(events, "L1"):
          H1L1_eff_snr = math.sqrt((events.H1.get_effective_snr())**2 + (events.L1.get_effective_snr())**2)
          return {
                "H1L1_eff_snr": H1L1_eff_snr
          }

        elif hasattr(events, "H2") and hasattr(events, "L1"):
          H2L1_eff_snr = math.sqrt((events.H2.get_effective_snr())**2 + (events.L1.get_effective_snr())**2)
          return {
                "H2L1_eff_snr": H2L1_eff_snr
          }


statistic = CoincInspiralUtils.coincStatistic(opts.statistic)

###############################################################################
# read in zero lag coinc triggers

zerolagTriggers = None
zerolagTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles(zerolagfiles, mangle_event_id = True)

# construct the zero lag coincs
zerolagCoincTriggers= \
CoincInspiralUtils.coincInspiralTable(zerolagTriggers, statistic)

slidesTriggers = None
slidesTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles(slidesfiles, mangle_event_id = True)

# construct the time slides coincs for single.
slidesCoincs= \
CoincInspiralUtils.coincInspiralTable(slidesTriggers, statistic)


# Construct the time slides for double in triple times.
if opts.coincs == "H1H2":
  H1H2_slides_Coincs = CoincInspiralUtils.coincInspiralTable(slidesTriggers, statistic).coincinclude(ifolist)
  
if opts.coincs == "H1L1": 
  H1L1_slides_Coincs = CoincInspiralUtils.coincInspiralTable(slidesTriggers, statistic).coinctype(ifolist)

if opts.coincs == "H2L1":
  H2L1_slides_Coincs = CoincInspiralUtils.coincInspiralTable(slidesTriggers, statistic).coinctype(ifolist)



# read in the injection files
injectionTriggers = None
injectionTriggers = SnglInspiralUtils.ReadSnglInspiralFromFiles(foundfiles, mangle_event_id = True)

# contruct the injection coincs for double in tripe times.

if opts.coincs == "H1H2":
  H1H2_Inj_Coincs = CoincInspiralUtils.coincInspiralTable(injectionTriggers, statistic).coincinclude(ifolist)

if opts.coincs == "H1L1":
  H1L1_Inj_Coincs = CoincInspiralUtils.coincInspiralTable(injectionTriggers, statistic).coinctype(ifolist)

if opts.coincs == "H2L1":
  H2L1_Inj_Coincs = CoincInspiralUtils.coincInspiralTable(injectionTriggers, statistic).coinctype(ifolist)

# Construct coincs for singles.

InjectionCoincs = CoincInspiralUtils.coincInspiralTable(injectionTriggers, statistic)

################################################################################
#
# Create a book-keeping object.  You initialize this thing by passing a
# series of key-word arguments.  The name of each argument becomes the name
# of one of the parameters whose distribution will be measured.  The value
# you set each argument to is a tuple to be used as the arguments of
# pylal.rate.BinnedArray instances created to track the "background" and
# "injections" distributions for that parameter.

###############################################################################

# Create a list of doubles and Triples
##############################################################################


#ifos = [opts.ifo_first, opts.ifo_second]
#ifos.sort()

#Ethinca_Injections = InjectionCoincs.coincinclude(ifos).getEThincaValues(ifos)     
#Ethinca_Slides = slidesCoincs.coincinclude(ifos).getEthincaValues(ifos)

#print Ethinca_Slides
###############################################################################


class CoincParamsDistributions(ligolw_burca_tailor.BurcaCoincParamsDistributions):
	binnings = {
        	"H1_eff_snr": rate.NDBins((rate.LinearBins(0.0, 50.0, 1000),)),
        	"H2_eff_snr": rate.NDBins((rate.LinearBins(0.0, 50.0, 1000),)),
        	"L1_eff_snr": rate.NDBins((rate.LinearBins(0.0, 50.0, 1000),)),
        	"H1H2_eff_snr": rate.NDBins((rate.LinearBins(0.0, 50.0, 1000),)),
        	"H1L1_eff_snr": rate.NDBins((rate.LinearBins(0.0, 50.0, 1000),)),
        	"H2L1_eff_snr": rate.NDBins((rate.LinearBins(0.0, 50.0, 1000),))
	}
	filters = {
        	"H1_eff_snr": rate.gaussian_window(21),
        	"H2_eff_snr": rate.gaussian_window(21),
        	"L1_eff_snr": rate.gaussian_window(21),
        	"H1H2_eff_snr": rate.gaussian_window(21),
        	"H1L1_eff_snr": rate.gaussian_window(21),
        	"H2L1_eff_snr": rate.gaussian_window(21)
	}
distributions = CoincParamsDistributions()

timeslide = 0


#############################################################################
# Constructions of effective snr arrays for Injections
############################################################################


x_inj_param = []
for coincs in InjectionCoincs:
#############################################################################
# For Single ifos
#############################################################################
  
  if opts.statistic == 'effective_snr' and opts.ifo == "H1" and hasattr(coincs, "H1"):
        distributions.add_injection(single_params_func(coincs.H1, timeslide))
        x_inj_param.append(coincs.H1.get_effective_snr())
 
  elif  opts.statistic == 'effective_snr' and opts.ifo == "H2" and hasattr(coincs, "H2"):
        distributions.add_injection(single_params_func(coincs.H2, timeslide))
        x_inj_param.append(coincs.H2.get_effective_snr())

  elif  opts.statistic == 'effective_snr' and opts.ifo == "L1" and hasattr(coincs, "L1"):
        distributions.add_injection(single_params_func(coincs.L1, timeslide))
        x_inj_param.append(coincs.L1.get_effective_snr())

##############################################################################
# For double in Double times
##############################################################################
  
  if opts.statistic == 'effective_snr' and opts.coincs == "H1H2" and hasattr(coincs, "H1") and hasattr(coincs, "H2"):
        distributions.add_injection(double_params_func(coincs, timeslide))
        x_inj_param.append(math.sqrt((coincs.H1.get_effective_snr())**2 + (coincs.H2.get_effective_snr())**2))
     
  elif opts.statistic == 'effective_snr' and opts.coincs == "H1L1" and hasattr(coincs, "H1") and hasattr(coincs, "L1"):
        distributions.add_injection(double_params_func(coincs, timeslide))
        x_inj_param.append(math.sqrt((coincs.H1.get_effective_snr())**2 + (coincs.L1.get_effective_snr())**2))

  elif opts.statistic == 'effective_snr' and opts.coincs == "H2L1" and hasattr(coincs, "H2") and hasattr(coincs, "L1"):     
        distributions.add_injection(double_params_func(coincs, timeslide))
        x_inj_param.append(math.sqrt((coincs.H2.get_effective_snr())**2 + (coincs.L1.get_effective_snr())**2))
          
X_inj_param = asarray(x_inj_param)
#############################################################################   
#Constructions of effective snr arrays for background
#############################################################################

x_back_param = []
for coincs in slidesCoincs:
  
  if opts.statistic == 'effective_snr' and opts.ifo == "H1" and hasattr(coincs, "H1"):
        distributions.add_background(single_params_func(coincs.H1, timeslide))
        x_back_param.append(coincs.H1.get_effective_snr())

  elif  opts.statistic == 'effective_snr' and opts.ifo == "H2" and hasattr(coincs, "H2"):
        distributions.add_background(single_params_func(coincs.H2, timeslide))
        x_back_param.append(coincs.H2.get_effective_snr())
  
  elif  opts.statistic == 'effective_snr' and opts.ifo == "L1" and hasattr(coincs, "L1"):
        distributions.add_background(single_params_func(coincs.L1, timeslide))
        x_back_param.append(coincs.L1.get_effective_snr())

#############################################################################
# Doubles in Triple times
#############################################################################

  if opts.statistic == 'effective_snr' and opts.coincs == "H1H2" and hasattr(coincs, "H1") and hasattr(coincs, "H2"):
        distributions.add_injection(double_params_func(coincs, timeslide))
        x_back_param.append(math.sqrt((coincs.H1.get_effective_snr())**2 + (coincs.H2.get_effective_snr())**2))

  elif opts.statistic == 'effective_snr' and opts.coincs == "H1L1" and hasattr(coincs, "H1") and hasattr(coincs, "L1"):
        distributions.add_injection(double_params_func(coincs, timeslide))
        x_back_param.append(math.sqrt((coincs.H1.get_effective_snr())**2 + (coincs.L1.get_effective_snr())**2))

  elif opts.statistic == 'effective_snr' and opts.coincs == "H2L1" and hasattr(coincs, "H2") and hasattr(coincs, "L1"):
        distributions.add_injection(double_params_func(coincs, timeslide))
        x_back_param.append(math.sqrt((coincs.H2.get_effective_snr())**2 + (coincs.L1.get_effective_snr())**2))


  
X_back_param = asarray(x_back_param)
#############################################################################
# Finish Smoothening of the Data using Gaussian Filter
#############################################################################
xmldoc = ligolw_burca_tailor.gen_likelihood_control(distributions)
utils.write_filename(xmldoc, "distributions.xml")   

distributions.finish()
#3###########################################################################
# Construction of Histrogram
#############################################################################

p = arange(0, 50.0, 1.0)
X_Inj_norm = hist(X_inj_param,p)[0]*1.0/max(hist(X_back_param,p)[0])
clf()

X_Back_norm = hist(X_back_param,p)[0]*1.0/max(hist(X_back_param,p)[0])
clf()

##########################################################################
# Reload X nd Y parameters
##########################################################################

# For Singles
##########################################################################

if opts.statistic == 'effective_snr' and opts.ifo == "H1":
  x_inj = distributions.injection_pdf["H1_eff_snr"].centres()[0]
  y_inj = distributions.injection_pdf["H1_eff_snr"].array

elif opts.statistic == 'effective_snr' and opts.ifo == "H2":
  x_inj = distributions.injection_pdf["H2_eff_snr"].centres()[0]
  y_inj = distributions.injection_pdf["H2_eff_snr"].array

elif  opts.statistic == 'effective_snr' and opts.ifo == "L1":
  x_inj = distributions.injection_pdf["L1_eff_snr"].centres()[0]
  y_inj = distributions.injection_pdf["L1_eff_snr"].array

#########################################################################

if opts.statistic == 'effective_snr' and opts.ifo == "H1":
  x_back = distributions.background_pdf["H1_eff_snr"].centres()[0]
  y_back = distributions.background_pdf["H1_eff_snr"].array
  print "y_back", y_back

elif opts.statistic == 'effective_snr' and opts.ifo == "H2":
  x_back = distributions.background_pdf["H2_eff_snr"].centres()[0]
  y_back = distributions.background_pdf["H2_eff_snr"].array
  

elif opts.statistic == 'effective_snr' and opts.ifo == "L1":
  x_back = distributions.background_pdf["L1_eff_snr"].centres()[0]
  y_back = distributions.background_pdf["L1_eff_snr"].array

########################################################################
# For doubles in Triple times.
########################################################################

if opts.statistic == 'effective_snr' and opts.coincs == "H1H2":
  x_inj = distributions.injection_pdf["H1H2_eff_snr"].centres()[0]
  y_inj = distributions.injection_pdf["H1H2_eff_snr"].array

elif opts.statistic == 'effective_snr' and opts.coincs == "H1L1":
  x_inj = distributions.injection_pdf["H1L1_eff_snr"].centres()[0]
  y_inj = distributions.injection_pdf["H1L1_eff_snr"].array

elif opts.statistic == 'effective_snr' and opts.coincs == "H2L1":
  x_inj = distributions.injection_pdf["H2L1_eff_snr"].centres()[0]
  y_inj = distributions.injection_pdf["H2L1_eff_snr"].array

########################################################################
# Adding time slides.
########################################################################

if opts.statistic == 'effective_snr' and opts.coincs == "H1H2":
  x_back = distributions.background_pdf["H1H2_eff_snr"].centres()[0]
  y_back = distributions.background_pdf["H1H2_eff_snr"].array

elif opts.statistic == 'effective_snr' and opts.coincs == "H1L1":
  x_back = distributions.background_pdf["H1L1_eff_snr"].centres()[0]
  y_back = distributions.background_pdf["H1L1_eff_snr"].array

  
elif opts.statistic == 'effective_snr' and opts.coincs == "H2L1":
  x_back = distributions.background_pdf["H2L1_eff_snr"].centres()[0]
  y_back = distributions.background_pdf["H2L1_eff_snr"].array

########################################################################
# Plottings routines
########################################################################

if (opts.ifo == "H1" or opts.ifo == "H2" or opts.ifo == "L1"):
  plot(x_inj, y_inj, "r-", x_back, y_back, "k-")
  title('Normalised Distribution of Background and Injections')
  legend(("Inj","Back"),loc=0)
  xlabel('Effective_snr' + " " + str(opts.ifo))
  ylabel('Counts')
  figure()
  bar(p, X_Inj_norm,color='r')
  hold(True)
  bar(p, X_Back_norm,color='k')
  legend(("Inj","Back"))
  title('Normalised histogram of Background and Injections')
  xlabel('Effective_snr' + " " + str(opts.ifo))
  ylabel('Counts')


if (opts.coincs == "H1H2" or opts.coincs == "H1L1" or opts.coincs == "H2L1"):
  plot(x_inj, y_inj, "r", x_back, y_back, "k-")
  title('Normalised Distribution of Background and Injections')
  legend(("Inj","Back"),loc=0)
  xlabel('Combined_snr' + " " + str(opts.coincs))
  ylabel('Counts')
  figure()
  bar(p, X_Inj_norm,color='r')
  hold(True)
  legend("Inj")
  bar(p, X_Back_norm,color='k')
  hold(True)
  legend("Back",loc = 0)
  title('Normalised histogram of Background and Injections')
  xlabel('Combined_snr' + " " + str(opts.coincs))
  ylabel('Counts')

  
########################################################################
show()
