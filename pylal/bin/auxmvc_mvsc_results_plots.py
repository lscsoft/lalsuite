#!/usr/bin/python
from optparse import *
import glob
import sys
import os
import matplotlib
matplotlib.use('Agg')
import pylab
import pdb
import numpy
from pylal import auxmvc_utils 
import bisect
import pickle

parser=OptionParser(usage="Generates summary plots for KW aux-triggers processed by StatPatternRecognition", version = "Kari Hodge")
parser.add_option("","--histograms", action="store_true", default=False, help="use if you want to produce histograms for each dimension")
parser.add_option("","--plot-rank-vs-significance", action="store_true", default=False, help="Make scatter plots of rank vs Aux KW significance")
parser.add_option("","--plot-rank-vs-dt", action="store_true", default=False, help="Make scatter plots of rank vs dt")
parser.add_option("","--tag", help="filenames will be ROC_tag.png and efficiency_deadtime_tag.txt")
parser.add_option("","--output-dir", help="directory where output files will be written to")
(opts,files)=parser.parse_args()

try: os.mkdir(opts.output_dir)
except: pass

print files
data = auxmvc_utils.ReadMVSCTriggers(files)        

clean_data = data[numpy.nonzero(data['i']==0)[0],:]
glitch_data = data[numpy.nonzero(data['i']==1)[0],:]

# The rank name is 'Bagger' for MVSC, 'glitch-rank' for ANN, and 'SVMRank' for SVM. The last column in *.dat files is ranking values, so the end component of first line in *.dat file is the rank name.
rank_name = data.dtype.names[-1]


all_ranks = numpy.concatenate((clean_data[rank_name],glitch_data[rank_name]))
all_ranks_sorted = numpy.sort(all_ranks)


FAP, TAP = auxmvc_utils.ROC(clean_data[rank_name], glitch_data[rank_name])
#FAP = false alarm percentage = number of random/clean times flagged as glitches
#TAP = true alarm percentage = number of glitches flagged as glitches

# Plot ROC curve
pylab.figure(1)
pylab.loglog(FAP,TAP, linewidth = 2.0)
pylab.hold(True)
x = numpy.arange(min(TAP), max(TAP) + (max(TAP) - min(TAP))/1000.0, (max(TAP) - min(TAP))/1000.0)
pylab.loglog(x,x, linestyle="dashed", linewidth = 2.0)
pylab.xlabel('False Alarm Probability')
pylab.ylabel('Efficiency')
pylab.xlim([0,1])
pylab.ylim([0,1])
pylab.savefig(opts.output_dir + "/"+'ROC_'+opts.tag+'.png')	
pylab.close()

# save ROC curve in a file
roc_file = open(opts.output_dir + "/"+"ROC_" + opts.tag + ".pickle", "w")
pickle.dump([FAP,TAP], roc_file)
roc_file.close()

#FAP is a list that is naturally sorted in reverse order (highest to lowest),
#we need to turn it into a regularly sorted list so that we can find the TAP for
#fiducial FAPs
FAP.sort()
edfile = open(opts.output_dir + "/"+'efficiency_deadtime_'+opts.tag+'.txt','w')
for threshold in [.01,.05,.1]:
	tmpindex=bisect.bisect_left(FAP,threshold)
	edfile.write("deadtime: "+str(FAP[tmpindex])+" efficiency: "+str(TAP[len(FAP)-tmpindex-1])+"\n")


# plot rank vs significance

if opts.plot_rank_vs_significance:

  variables = data.dtype.names
  #sig_variables = [s for s in variables if s.endswith("_sig")]
  sig_variables = [s for s in variables if s.endswith("_signif")]
  #form arrays with significances for clean data sets
  clean_sig_rss = []
  clean_sig_max = []
  clean_rank = []
  for i in range(len(clean_data)):
    sig_rss = 0
    sig_max = 0
    for var in sig_variables:
      sig_rss += clean_data[var][i]**2
      if clean_data[var][i] > sig_max:
        sig_max = clean_data[var][i]
    clean_sig_rss.append(sig_rss**0.5)
    clean_sig_max.append(sig_max)
    clean_rank.append(clean_data[rank_name][i])
    

  #form arrays with significances for glitch data sets
  glitch_sig_rss = []
  glitch_sig_max = []
  glitch_rank = []
  for i in range(len(glitch_data)):
    sig_rss = 0
    sig_max = 0
    for var in sig_variables:
      sig_rss += glitch_data[var][i]**2
      if glitch_data[var][i] > sig_max:
        sig_max = glitch_data[var][i]
    glitch_sig_rss.append(sig_rss**0.5)
    glitch_sig_max.append(sig_max)
    glitch_rank.append(glitch_data[rank_name][i])
    

  pylab.figure(1)
  pylab.loglog(clean_sig_max,clean_rank, "kx", label="clean")
  pylab.hold(True)
  pylab.loglog(glitch_sig_max,glitch_rank, "r+", label="glitch")
  pylab.hold(True)
  pylab.xlabel('Max Aux KW significance')
  pylab.ylabel('MVSC rank')
  pylab.savefig(opts.output_dir + "/"+'rank_vs_sig_max_'+opts.tag+'.png')	
  pylab.close()


  pylab.figure(1)
  pylab.loglog(clean_sig_rss,clean_rank, "kx", label="clean")
  pylab.hold(True)
  pylab.loglog(glitch_sig_rss,glitch_rank, "r+", label="glitch")
  pylab.hold(True)
  pylab.xlabel('rss Aux KW significance')
  pylab.ylabel('MVSC rank')
  pylab.savefig(opts.output_dir + "/"+'rank_vs_sig_rss_'+opts.tag+'.png')	
  pylab.close()


# plot rank vs dt

if opts.plot_rank_vs_dt:

  variables = data.dtype.names
  sig_variables = [s for s in variables if s.endswith("_signif")]
  dt_variables = [s for s in variables if s.endswith("_dt")]
  
  #form arrays with dt for clean data sets
  clean_dt_min = []
  clean_dt_sig_max = []
  clean_rank = []
  for i in range(len(clean_data)):
    dt_min = numpy.inf
    sig_max = 0
    for (dt_var,sig_var) in zip(dt_variables, sig_variables):
      if (clean_data[sig_var][i] > 0) and (abs(clean_data[dt_var][i]) < dt_min):
        dt_min = abs(clean_data[dt_var][i])
      if clean_data[sig_var][i] > sig_max:
        dt_sig_max = abs(clean_data[dt_var][i])
      
      
    clean_dt_min.append(dt_min)
    clean_dt_sig_max.append(dt_sig_max)
    clean_rank.append(clean_data[rank_name][i])



  #form arrays with dt for glitch data sets
  glitch_dt_min = []
  glitch_dt_sig_max = []
  glitch_rank = []
  for i in range(len(glitch_data)):
    dt_min = numpy.inf
    sig_max = 0
    for (dt_var,sig_var) in zip(dt_variables, sig_variables):
      if (glitch_data[sig_var][i] > 0) and (abs(glitch_data[dt_var][i]) < dt_min):
        dt_min = abs(glitch_data[dt_var][i])
      if glitch_data[sig_var][i] > sig_max:
        dt_sig_max = abs(glitch_data[dt_var][i])
      
      
    glitch_dt_min.append(dt_min)
    glitch_dt_sig_max.append(dt_sig_max)
    glitch_rank.append(glitch_data[rank_name][i])
    
    

  pylab.figure(1)
  pylab.loglog(clean_dt_min,clean_rank, "kx", label="clean")
  pylab.hold(True)
  pylab.loglog(glitch_dt_min,glitch_rank, "r+", label="glitch")
  pylab.hold(True)
  pylab.xlabel('min dt, seconds')
  pylab.ylabel('MVSC rank')
  pylab.savefig(opts.output_dir + "/"+'rank_vs_dt_min_'+opts.tag+'.png')	
  pylab.close()

  pylab.figure(1)
  pylab.loglog(clean_dt_sig_max, clean_rank, "kx", label="clean")
  pylab.hold(True)
  pylab.loglog(glitch_dt_sig_max, glitch_rank, "r+", label="glitch")
  pylab.hold(True)
  pylab.xlabel('dt for trigger with highest sig, seconds')
  pylab.ylabel('MVSC rank')
  pylab.savefig(opts.output_dir + "/"+'rank_vs_dt_sig_max_'+opts.tag+'.png')	
  pylab.close()


if opts.histograms:
  variables = data.dtype.names
  fig_num = 1000
  for i,var in enumerate(variables):
    fig_num += 1
    pylab.figure(fig_num)
    print var
    pylab.hist(clean_data[var],100)
    pylab.savefig(opts.output_dir + "/"+"hist_twosided"+var+"_cleantimes") 
    pylab.close()
    fig_num += 1
    pylab.figure(fig_num)
    pylab.hist(glitch_data[var],bins=100)
    pylab.savefig(opts.output_dir + "/"+"hist_twosided"+var+"_glitches")
    pylab.close() 
