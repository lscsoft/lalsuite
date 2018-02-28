#!/usr/bin/python
# Reed Essick (reed.essick@ligo.org)

#written to load pickle files and generate plots accordingly

__author__ = "Reed Essick"
__version__ = "$Revision$"[11:-2]
__date__ = "$Date$"[7:-2]

__prog__="auxmvc_comparison_plots_from_pickle_files.py"
__Id__ = "$Id$"



from optparse import *
import glob
import sys
import os
import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.axes_grid import ImageGrid
import pylab
import pdb
import numpy
from pylal import auxmvc_utils 
from pylal import git_version
import bisect
import pickle
from pylal import InspiralUtils
import math

usage= """Written to load and manipulate output from different classifiers."""

parser=OptionParser(usage=usage, version = git_version.verbose_msg)

parser.add_option("","--fap-threshold", default=0.1,type="float", help="False Alarm Probability which is adapted to veto")
parser.add_option("-P","--output-path",action="store",type="string",default="",  metavar="PATH", help="path where the figures would be stored")	  
parser.add_option("-O","--enable-output",action="store_true", default="True",  metavar="OUTPUT", help="enable the generation of the html and cache documents")	 	  
parser.add_option("-u","--user-tag",action="store",type="string", default=None,metavar=" USERTAG", help="a user tag for the output filenames" )
parser.add_option("", "--figure-resolution",action="store",type="int", default=50, help="dpi of the thumbnails (50 by default)")
parser.add_option("", "--html-for-cbcweb",action="store", default=False, metavar = "CVS DIRECTORY", help="publish the html "\
      "output in a format that can be directly published on the cbc webpage "\
      "or in CVS. This only works IF --enable-output is also specified. The "\
      "argument should be the cvs directory where the html file will be placed "\
      "Example: --html-for-cbcweb protected/projects/s5/yourprojectdir")
parser.add_option("","--verbose", action="store_true", default=False, help="print information" )

# plotting options
parser.add_option("","--DQ-ROC", type='string', default=False, help='plots DQ flags on ROC plots. DQ-ROC can be either S4 or S6') 
parser.add_option("","--diagnostic-plots", action="store_true", default=False, help="generates a large number of diagnostic plots meant to help determine the relative performance and correlation of classifiers")


# plotting pickle file locations
parser.add_option("","--ROC", default=False, type="string", help="Provide the path for the ROC pickle file")
parser.add_option("","--bit-word", default=False, type="string", help="Provide the path for the bit-word pickle file")
parser.add_option("","--hist-ranks", default=False, type="string", help="Provide the path for the hist-ranks pickle file")
parser.add_option("","--cum-hist-signif", default=False, type="string", help="Provide the path for the cum-hist-signif pickle file")

# use glob.glob('filename') to extract files from opts.pickle_files

(opts,args)=parser.parse_args()

try: os.mkdir(opts.output_path)
except: pass


################   PLOTS   #############################################################

# Initializing the html output
InspiralUtils.message(opts, "Initialisation...")
opts = InspiralUtils.initialise(opts, __prog__, __version__)
fnameList = []
tagList = []
fig_num = 0
comments = ""

###########################################################################################
if opts.verbose:
  print 'Generating Plots...'

colorDIC = {'ann':'b', 'mvsc':'g', 'svm':'r', 'ovl':'c', 'combined':'m'}
labelDIC = {'ann':'ANN', 'mvsc':'MVSC', 'svm':'SVM', 'ovl':'OVL', 'combined':'MVC$_{\mathrm{max}}$'}

matplotlib.rc('text', usetex=True)

faircoin = numpy.linspace(10**-6,10**0.5,100)


#############################################################################################
#
### ROC curves
#
# generates the standard ROC curves using values for eff and fap computed within this script.
#
#############################################################################################

if opts.ROC:

  if opts.verbose:
    print '  ROC curves'  

  # open and load pickle file for ROC data
  pfile = open(opts.ROC) # the order in which data was dumped for ROC pickle files is:
                                    #   classifiers
                                    #   cls1_fap
                                    #   cls1_eff
                                    #   cls2_fap
                                    #     ...
                                    #   clsN_fap
                                    #   clsN_eff
                                    # where N = len(classifiers)

  classifiers = pickle.load(pfile) # labels the data in the order it was dumped

  # generate figure
  fig_num += 1
  fig = pylab.figure(fig_num)
  for cls in classifiers:
    fap = pickle.load(pfile)
    eff = pickle.load(pfile)
    pylab.plot(fap, eff, label = labelDIC[cls[0]], color=colorDIC[cls[0]])#, linewidth = 2)
  pylab.plot(faircoin, faircoin, '--k')
  pylab.xlabel('False Alarm Probability')
  pylab.ylabel('Efficiency')
  pylab.xscale('log')
  pylab.yscale('log')
  pylab.xlim(10**-5, 10**0)
  pylab.ylim(10**-3.5, 10**0)
  pylab.title('Fig. '+str(fig_num)+': ROC Curves from total\_ranked\_data[]')
  l1 = pylab.legend(loc = 'lower right')

  pfile.close()

  if opts.diagnostic_plots:
    # save log-log figure
    name = '_ROC_log-log'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)

    # save lin-lin figure
    fig_num +=1
    pylab.title('Fig. '+str(fig_num)+': ROC Curves from total\_ranked\_data[]')
    pylab.legend(loc = 'lower right')
    pylab.yscale('linear')
    pylab.xscale('linear')
    pylab.xlim(0, 1)
    pylab.ylim(0, 1)
    name = '_ROC_lin-lin'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)

    fig_num += 1

  # check to see which DQ_ROC points to plot
  if opts.DQ_ROC == 'S4':
    cbcS4DQfap = [2355./98147, 7672./98147, 19179./98147]
    cbcS4DQeff = [2357./16205, 6710./16205, 9445./16205]
    dq2,   = pylab.plot(cbcS4DQfap[0], cbcS4DQeff[0], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = '^', markersize = 6, label = 'CBC DQ\,2')
    dq23,  = pylab.plot(cbcS4DQfap[1], cbcS4DQeff[1], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 's', markersize = 6, label = 'CBC DQ\,2+3')
    dq234, = pylab.plot(cbcS4DQfap[2], cbcS4DQeff[2], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 'd', markersize = 6, label = 'CBC DQ\,2+3+4')
    pylab.legend(loc = 'upper left')
  
  if opts.DQ_ROC == 'S6':
    burstS6DQfap = [485./99869, 4028./99869, 4503./99869]
    burstS6DQeff = [547./2833, 721./2833, 748./2833]
    dq2,   = pylab.plot(burstS6DQfap[0], burstS6DQeff[0], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = '^', markersize = 6, label = 'Burst DQ\,2')
    dq23,  = pylab.plot(burstS6DQfap[1], burstS6DQeff[1], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 's', markersize = 6, label = 'Burst DQ\,2+3')
    dq234, = pylab.plot(burstS6DQfap[2], burstS6DQeff[2], linestyle = 'None', markerfacecolor = 'None', markeredgecolor = 'k', marker = 'd', markersize = 6, label = 'Burst DQ\,2+3+4')
    pylab.legend(loc = 'upper left')

  # save lin-log figure  
  pylab.title(opts.user_tag + '\nROC Curves')
  pylab.legend(loc = 'upper left', ncol=2)
  pylab.xscale('log')
  pylab.yscale('linear')
  pylab.xlim(10**-5, 10**0)
  pylab.ylim(0,1)
  	
  name = '_ROC_lin-log'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()


#############################################################################################
#
### Bit-word histograms
#
# generate histograms representing the fractions of glitches removed by different combinations of classifiers
# 
#
#############################################################################################

if opts.bit_word:

  if opts.verbose:
    print '  bit-word histograms'

  # open pickle file
  pfile = open(opts.bit_word) # pickle file was dumped in this order:
                              #  sigthr
                              #  FAPthr
                              #  clas (labels the order of the bit-word)
                              #  fbw (all found glitches)
                              #  fbwBelow
                              #  fbwAbove

  sigthr = pickle.load(pfile)
  FAPthr = pickle.load(pfile)

  ### we define a 5 bit word to be the concatenation of binary flags corresponding to whether or not a classifier removes a glitch at a set FAP
  # 0: classifier did not remove glitch
  # 1: classifier removed glitch

  ### this next bit is fragile, and only handles the special case of svm being absent. you should extend this to a more general case

  clas = pickle.load(pfile)

  fbw = pickle.load(pfile)
  fbwBelow = pickle.load(pfile)
  fbwAbove = pickle.load(pfile)
  
  pfile.close()
    
  if opts.diagnostic_plots:
    ### Add regular histograms 
    fig_num += 1
    fig = matplotlib.pyplot.figure()
    pylab.hist(fbw, bins = numpy.linspace(-0.5, 2**len(clas)-0.5, num = 2**len(clas) +1, endpoint=True), weights = numpy.ones(len(fbw))/len(fbw), histtype='step', log = True, label = 'all glitches')#, linewidth = 2)
    pylab.hist(fbwBelow, bins = numpy.linspace(-0.5, 2**len(clas)-0.5, num = 2**len(clas) +1, endpoint = True), weights = numpy.ones(len(fbwBelow))/len(fbwBelow), histtype = 'step', log = True, label = 'signif $\leq$ ' + str(sigthr))#,  linewidth = 2)
    (numA, bins, pathces) = pylab.hist(fbwAbove, bins = numpy.linspace(-0.5, 2**len(clas)-0.5, num = 2**len(clas) +1, endpoint = True), weights = numpy.ones(len(fbwAbove))/len(fbwAbove), histtype = 'step', log = True, label = 'signif $\geq$ ' + str(sigthr))#, linewidth = 2)
    pylab.ylabel('Fraction of Glitches')
    pylab.xlim(-0.5, 2**len(clas)-0.5)
    pylab.ylim(ymax = 1.0)
    pylab.title('Fig. '+str(fig_num)+': Histogram over Classifier Redundancy at FAP = '+str(FAPthr))
    pylab.legend(loc = 'upper center')
    pylab.grid(True, which = 'major')
    pylab.grid(True, which = 'minor')
   
    #convert decimals back to binaries and label the ticks accordingly
    tick_locs = numpy.linspace(0,2**len(clas)-1,num=2**len(clas))
    tick_labels = []
    for loc in tick_locs:
      s = ''
      l = range(len(clas) - 1)
      for ind in l[::-1]:
        if int(loc)/2**(ind+1) == 0:
          s += '0'
        else:
          s += '1'
          loc = int(loc)-2**(ind+1)
      if int(loc) == 0:
        s += '0'
      else:
        s += '1'
      tick_labels += [s]
    pylab.xticks(tick_locs, tick_labels)
    for label in fig.axes[0].xaxis.get_ticklabels():
      label.set_rotation(45)
      leg = "5 bit word is ordered in the following way:  "
    for cls in clas:
       leg += labelDIC[cls[0]]+' '
    pylab.figtext(0.5, 0.98, leg, ha = 'center', va = 'top')
    
    # print the fractions on top of the steps in the histograms
    for indn in range(len(numA)):
      loc = tick_locs[indn]
      if numA[indn] > 0:
        n = int(numA[indn]*10**3)/10.0
        pylab.text(loc,numA[indn],str(n)+'\,\%',ha='center',va='bottom')
  
    #adding to html page
    strclas = ''
    for cls in clas:
      strclas += cls[0] + '_'
    name = '_scatter_5bit-words_'+strclas+'_FAPthr_'+str(FAPthr)
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()
  
  ### add histograms over only those glitches that were removed
  fig_num += 1
  fig = matplotlib.pyplot.figure()
  fbwFound = []
  for fbwi in fbw:
    if fbwi > 0:
      fbwFound += [fbwi]
  fbwAboveFound = []
  for fbwA in fbwAbove:
    if fbwA > 0:
      fbwAboveFound += [fbwA]
  fbwBelowFound = []
  for fbwB in fbwBelow:
    if fbwB > 0:
      fbwBelowFound += [fbwB]

  pylab.hist(fbwFound, bins = numpy.linspace(0.5, 2**len(clas)-0.5, num = 2**len(clas), endpoint=True), weights = numpy.ones(len(fbwFound))/len(fbwFound), histtype='step', log = True, label = 'all glitches')#, linewidth = 2)
  pylab.hist(fbwBelowFound, bins = numpy.linspace(0.5, 2**len(clas)-0.5, num = 2**len(clas), endpoint = True), weights = numpy.ones(len(fbwBelowFound))/len(fbwBelowFound), histtype = 'step', log = True, label = 'signif $\leq$ ' + str(sigthr))#, linewidth = 2)
  (numA, bins, patches) = pylab.hist(fbwAboveFound, bins = numpy.linspace(0.5, 2**len(clas)-0.5, num = 2**len(clas), endpoint = True), weights = numpy.ones(len(fbwAboveFound))/len(fbwAboveFound), histtype = 'step', log = True, label = 'signif $\geq$ ' + str(sigthr))#, linewidth = 2)
  pylab.ylabel('Fraction of Glitches Found')
  pylab.xlim(0.5, 2**len(clas)-0.5)
  pylab.ylim(ymax = 1.0, ymin = 0.0)
  pylab.yscale('linear')
  pylab.title(opts.user_tag + '\nHistogram over Classifier Redundancy at FAP = '+str(FAPthr))
  pylab.legend(loc = 'upper right')
  pylab.grid(True, which = 'major')
  pylab.grid(True, which = 'minor')
  
  #convert decimals back to binaries and label the ticks accordingly
  tick_locs = numpy.linspace(1,2**len(clas)-1,num=2**len(clas)-1)
  tick_labels = []
  for loc in tick_locs:
    s = ''
    l = range(len(clas) - 1)
    for ind in l[::-1]:
      if int(loc)/2**(ind+1) == 0:
        s += '0'
      else:
        s += '1'
        loc = int(loc)-2**(ind+1)
    if int(loc) == 0:
      s += '0'
    else:
      s += '1'
    tick_labels += [s]
  pylab.xticks(tick_locs, tick_labels)
  for label in fig.axes[0].xaxis.get_ticklabels():
    label.set_rotation(45)
   
  # print the fractions on top of the steps in the histograms
  for indn in range(len(numA)):
    loc = tick_locs[indn]
    if numA[indn] > 0:
      n = int(numA[indn]*10**3)/10.0
      pylab.text(loc,numA[indn],str(n)+'\,\%',ha='center',va='bottom')
  
  leg = "Bit-word ordering:\n("
  for tmp in range(len(clas)-1):
    cls = clas[tmp]
    leg += labelDIC[cls[0]]+',\ '
  leg += labelDIC[clas[len(clas)-1][0]]+')'
  pylab.figtext(0.3, 0.8, leg, ha = 'center', va = 'center')

  #adding to html page
  strclas = ''
  for cls in clas:
    strclas += cls[0] + '_'
  name = '_scatter_5bit-words_FOUND_'+strclas+'_FAPthr_'+str(FAPthr)
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()

  
#############################################################################################
#
### Overlay of histograms of glitches and histograms of cleans over cls_rank
#
#############################################################################################

if opts.hist_ranks:

  if opts.verbose:
    print '  Overlay of histogras of glitches and histograms of cleans over Classifier rank'

 # open pickle file
  pfile = open(opts.hist_ranks) # pickle file was dumped in this order:
                                #  cls (cls[0] labels which classifier was used)
                                #  cleans_ranks
                                #  counts (contains lists of ranks for glitches with different values of DARM signif 
                                #    (binned with partitions according to sigthrs)
                                #  labels (the labels corresponding to counts)
                                #  sigthrs
                                #  rankthr
                                #  fapthr

  cls = pickle.load(pfile)
  cleans = pickle.load(pfile)
  counts = pickle.load(pfile)
  labels = pickle.load(pfile)
  sigthrs = pickle.load(pfile)
  rankthr = pickle.load(pfile)
  fapthr = pickle.load(pfile)

  pfile.close()

  # BARSTACKED histograms
  fig_num += 1
  pylab.figure(fig_num)
  pylab.hold(True)

  pylab.hist(counts, 100, histtype = 'barstacked', label = labels)
  pylab.hist(cleans, 100, histtype='step', label = 'all cleans', log = True) #weights = numpy.ones(len(cleans[cls[0]+'_rank']))/len(cleans[cls[0]+'_rank']) )

  pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
  pylab.xlabel(labelDIC[cls[0]]+ '\_rank')
  pylab.ylabel('Number of Glitches')
  pylab.legend(loc = 'upper center')

  # plot a line corresponding to opts.fap_threshold
  lims = matplotlib.pyplot.axis()
  fapthrLINE = pylab.plot([rankthr, rankthr], [lims[2], lims[3]], color = 'k', linewidth = 2)
  fapthrTEXT = pylab.text(rankthr, 0.5*(lims[3]-lims[2])+lims[2], 'FAP = '+str(fapthr)+' \n'+cls[0]+'\_rank = ' + str(int(rankthr*10**4)/10**4.0) +' ', ha = 'right', va = 'center')
  matplotlib.pyplot.axis(lims)

  #adding to the html page
  name = '_hist_glitches_and_cleans_' + cls[0] + '_rank_LOG_barstacked'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)

  #adding to the html page
  pylab.yscale('linear')
  fig_num += 1
  pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
  name = '_hist_glitches_and_cleans' + cls[0] + '_rank_LINEAR_barstacked'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()

  # STEP histograms
  fig_num += 1
  pylab.figure(fig_num)
  pylab.hold(True)

  # histogram using only a subset of glitches with sufficiently large 'signif'
  g_rem = []
  for sigthrIND in range(len(sigthrs)): #[10, 15, 25, 50]:
    sigthr = sigthrs[(len(sigthrs) - 1) - sigthrIND]
    g_rem = g_rem + [i for i in counts[(len(counts)-1) - sigthrIND]]
    pylab.hist(g_rem, 100, histtype = 'step', weights = numpy.ones(len(g_rem))/len(g_rem), label = 'signif $\geq$ ' + repr(sigthr), log=True)

  # histogram over all glichtes
  g_rem = g_rem + [i for i in counts[0]]
  pylab.hist(g_rem, 100, histtype = 'step', weights = numpy.ones(len(g_rem))/len(g_rem), label = 'all glitches', log = True)

  pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
  pylab.xlabel(labelDIC[cls[0]]+ '\_rank')
  pylab.ylabel('Fraction of Glitches')

  # histogram using all the cleans
  pylab.hist(cleans, 100, histtype='step', weights = numpy.ones(len(cleans))/len(cleans), label = 'all cleans', log = True)

  pylab.legend(loc = 'upper center')

  # plot a line corresponding to opts.fap_threshold
  lims = matplotlib.pyplot.axis()
  fapthrLINE = pylab.plot([rankthr, rankthr], [lims[2], lims[3]], color = 'k', linewidth = 2)
  fapthrTEXT = pylab.text(rankthr, 0.5*(lims[3]-lims[2])+lims[2], 'FAP = '+str(fapthr)+' \n'+cls[0]+'\_rank = ' + str(int(rankthr*10**4)/10**4.0) +' ', ha = 'right', va = 'center')
  matplotlib.pyplot.axis(lims)

  #adding to the html page
  name = '_hist_glitches_and_cleans_' + cls[0] + '_rank_LOG'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)

  #adding to the html page
  pylab.yscale('linear')
  fig_num += 1
  pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
  name = '_hist_glitches_and_cleans_' + cls[0] + '_rank_LINEAR'
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()

"""
#############################################################################################
#
### Histograms of glitches over classifiers' ranks
#
# pretty self explanitory. We generate separate plots for just glitches, just cleans, and both glitches and cleans.
#
#############################################################################################

if opts.hist_ranks:

  if opts.verbose:
    print '  histograms of glitches over Classifier rank'

  for cls in classifiers:
    fig_num += 1
    pylab.figure(fig_num)
    pylab.hold(True)

    # histogram using all the glitches
    pylab.hist(glitches[cls[0] + '_rank'], 100, histtype='step', label = 'all glitches')

    # histogram using only a subset of glitches with sufficiently large 'signif'
    for sigthr in [10, 15, 25, 50]:
      glitches_removed = glitches[numpy.nonzero(glitches['signif'] >= sigthr)[0],:]
      if cls[0] == 'ovl':
        pylab.hist(glitches_removed[cls[0]+'_rank'], 100, histtype = 'step', label = 'signif $\geq$ ' + repr(sigthr), log=True)
      else:
        pylab.hist(glitches_removed[cls[0] + '_rank'], 100, histtype = 'step', label = 'signif $\geq$ ' + repr(sigthr), log=True)
    pylab.title('Fig. '+str(fig_num)+': Histogram for Glitches Based on ' + labelDIC[cls[0]] + '\_rank')
    pylab.xlabel(labelDIC[cls[0]]+ '\_rank')
    pylab.ylabel('Number of Glitches')
    pylab.legend(loc = 'upper center')

    #adding to the html page
    name = '_hist_glitches_' + cls[0] + '_rank_LOG'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)

    #adding to the html page
    pylab.yscale('linear')
    name = '_hist_glitches_' + cls[0] + '_rank_LINEAR'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()

#############################################################################################
#
### Histogras of clean samples over classifiers' ranks
#
#############################################################################################

if opts.hist_ranks:

  if opts.verbose:
    print '  histograms of clean samples over Classifier rank'
  for cls in classifiers:
    fig_num +=1
    pylab.figure(fig_num)
    pylab.hold(True)

    #histogram using all clean samples
    pylab.hist(cleans[cls[0]+'_rank'], 100, histtype='step',label = 'all', log = True)
    pylab.title('Fig. '+str(fig_num)+': Histogram for Clean Samples Based on '+labelDIC[cls[0]]+'\_rank')
    pylab.xlabel(labelDIC[cls[0]]+'\_rank')
    pylab.ylabel('Number of Samples')

    #adding to html page
    name = '_hist_cleans'+cls[0]+'_rank-LOG'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)

    #adding to html page
    pylab.yscale('linear')
    name = '_hist_cleans'+cls[0]+'_rank-LINEAR'
    fname = InspiralUtils.set_figure_name(opts, name)
    fname_thumb = InspiralUtils.savefig_pylal(filename = fname, doThumb = True, dpi_thumb=opts.figure_resolution)
    fnameList.append(fname)
    tagList.append(name)
    pylab.close()
"""


#############################################################################################
#
### Cumulative histograms of DARM significance for glitches before and after vetoing
#
#############################################################################################	

if opts.cum_hist_signif:

  pfile = open(opts.cum_hist_signif)
  FAPThr = pickle.load(pfile)
  classifiers = pickle.load(pfile)
  
  # histograms for SNR
  if opts.verbose:
    print '  cumulative histograms over signif'

  fig_num += 1
  pylab.figure(fig_num)
  g = pickle.load(pfile)
  pylab.hist(g, 400, histtype='step', cumulative =-1, label='before vetoing', color='k')
  for cls in classifiers:
    g_rem = pickle.load(pfile)
    pylab.hist(g_rem, 400, histtype='step', cumulative =-1, label=labelDIC[cls[0]], color=colorDIC[cls[0]]) 
  pylab.title(opts.user_tag + '\nCumulative Histogram of Glitch DARM Significance at FAP = '+str(FAPThr))
  pylab.xlabel('Significance')
  pylab.ylabel('Number of Glitches')
  pylab.xscale('log')
  pylab.yscale('log')
  pylab.xlim(xmin=min(g), xmax=4*10**2)
  pylab.legend()

  pfile.close()

  # adding to html page
  name = '_cumul_hist_signif_fap'+str(FAPThr)
  fname = InspiralUtils.set_figure_name(opts, name)
  fname_thumb = InspiralUtils.savefig_pylal(filename=fname, doThumb=True, dpi_thumb=opts.figure_resolution)
  fnameList.append(fname)
  tagList.append(name)
  pylab.close()


##############################################################################################################
### this next bit will help to print the options nicely on the html page
opts_dict = vars(opts)

html_filename = InspiralUtils.write_html_output(opts, [" --"+key.replace('_','-')+"="+str(opts_dict[key]) for key in opts_dict.keys()], fnameList, tagList, comment=comments)
InspiralUtils.write_cache_output(opts, html_filename, fnameList)

if opts.html_for_cbcweb:
  html_filename_publish = InspiralUtils.wrifacte_html_output(opts, args, fnameList, tagList, cbcweb=True)

##############################################################################################################
if opts.verbose:
  print 'Done.'

