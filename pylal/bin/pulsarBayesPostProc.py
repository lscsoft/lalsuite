#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       pulsarBayesPostProc.py
#
#       Copyright 2012
#       Matthew Pitkin <matthew.pitkin@ligo.org>,
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.
#

# This is based heavily on cbcBayesPostProc.py written by:
#       Benjamin Aylott <benjamin.aylott@ligo.org>,
#       Benjamin Farr <bfarr@u.northwestern.edu>,
#       Will M. Farr <will.farr@ligo.org>,
#       John Veitch <john.veitch@ligo.org>

#standard library imports
import sys
import os

from math import ceil,floor
import cPickle as pickle

from time import strftime

#related third party imports
import numpy as np

import matplotlib
matplotlib.use("Agg")

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version

from lalapps import pulsarpputils as pppu

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date


def pulsarBayesPostProc( outdir, data,
                         upperlimit, histbins,
                         #nested sampling options
                         priorfile, parfile, ns_Nlive,
                         # analysis data
                         Bkdata=None, ifos=None,
                         # pulsar information
                         corfile=None #,
                         #Turn on 2D kdes
                         #twodkdeplots=False,
                         #Turn on R convergence tests
                         #RconvergenceTests=False
                       ):
  # check data is input
  if data is None:
    raise RuntimeError('You must specify an input data file')
  
  # check output directory is input
  if outdir is None:
    raise RuntimeError("You must specify an output directory.")
  
  # create output dir if it doesn't exist
  if not os.path.isdir(outdir):
    os.makedirs(outdir)
  
  # list of posteriors for each ifo
  poslist = []
  
  # loop over ifos and generate posteriors
  for idx, ifo in enumerate(ifos):
    # nested sampling data files for IFO idx
    idat = data[idx]
    
    # if input format is xml
    votfile = None
    if '.xml' in data[0]:
      peparser = bppu.PEOutputParser('xml')
      commonResultsObj = peparser.parse(idat[0])
      thefile = open(idat[0],'r')
      votfile = thefile.read()
    else: # input as standard nested sampling input
      peparser = bppu.PEOutputParser('ns')
      commonResultsObj = peparser.parse(idat, Nlive=ns_Nlive, xflag=False)

    # Create an instance of the posterior class using the posterior values
    # loaded from the file.
    pos = bppu.Posterior( commonResultsObj, SimInspiralTableEntry=None, \
                          votfile=votfile )
   
    # append to list
    poslist.append(pos)
  
  # read in par file if given
  if parfile is not None:
    par = pppu.psr_par(parfile)
  
  # convert phi0' and psi' to phi0 and psi
  if 'phi0prime' in pos.names and 'psiprime' in pos.names:
    phi0samps, psisamps = pppu.phipsiconvert( pos['phi0prime'].samples, \
                                              pos['psiprime'].samples )
    # append phi0 and psi to posterior
    phi0_pos = bppu.OneDPosterior('phi0', phi0samps)
    pos.append(phi0_pos)
    
    psi_pos = bppu.OneDPosterior('psi', phi0samps)
    pos.append(psi_pos)
   
  # get parameters to plot from prior file and those with "Errors" in the par
  # file
  plotpars = []
  pri = pppu.psr_prior(priorfile)
  for pars in pppu.float_keys:
    # prior file values have preference over par file values with errors in the
    # nested sampling code
    pars_err = pars+'_ERR'
    
    if pri[pars] is not None:
      plotpars.append(pars)
    elif par[pars_err] is not None:
      plotpars.append(pars)
  
  # if h0 exists plot the pdf
  if 'H0' in plotpars:
    # set bounds of h0
    bounds = [0, float("inf")]
    
    plotFig, ulvals = pppu.plot_posterior_hist( poslist, 'H0'.lower(), ifos, 
                                                bounds, histbins, upperlimit )
    
    # output plots for each IFO
    for i, ifo in enumerate(ifos):
      figname = 'H0'.lower()+'_'+ifo+'.png'
      oneDplotPath = os.path.join(outdir, figname)
      plotFig[i].savefig(oneDplotPath)
      
  # if the heterodyned data files have been included create time series plots
  # of |B_k| and a amplitude spectral density spectrogram of the data
  if Bkdata is not None:
    Bkfigs, ASDfigs = pppu.plot_Bks_ASDs( Bkdata, ifos )
    
    # output plots
    for i, ifo in enumerate(ifos):
      figname = 'Bk'+'_'+ifo+'.png'
      Bkplotpath = os.path.join(outdir, figname)
      Bkfigs[i].savefig(Bkplotpath)
      
      figname = 'ASD'+'_'+ifo+'.png'
      ASDplotpath = os.path.join(outdir, figname)
      ASDfigs[i].savefig(ASDplotpath)
  
    
# main function
if __name__=='__main__':
  from optparse import OptionParser
  
  description = \
"""This script is for creating a results output page and plots for the known
   pulsar analysis. It uses inputs from the nested sampling code
   lalapps_pulsar_parameter_estimation_nested."""

  epilog = " An example of usage for a case when three nested sampling runs \
have been performed for two interferometers (H1 and L1): \
"+os.path.basename(sys.argv[0])+" --ifo H1 --data \
nest1_H1.txt,nest2_H1.txt,nest3_H1.txt --ifo L1 --data \
nest1_L1.txt,nest2_L1.txt,nest3_L1.txt --parfile J0534-2200.par --Nlive 1000 \
--priorfile priors.txt --histbins 50 --outpath \
/home/me/public_html/J0534-2200"
  usage = "Usage: %prog [options]"
  
  parser = OptionParser( usage = usage, description = description,
                         version = __version__, epilog = epilog )
  
  parser.add_option("-o", "--outpath", dest="outpath", help="The path for "
                    "the analysis output", metavar="DIR")
  
  """
   data will be read in in the following format:
     multiple nested sampling files from each IFO would be entered as follows:
     --ifo H1 --data nest1_H1.txt,nest2_H1.txt,nest3_H1.txt
     --ifo L1 --data nest1_L1.txt,nest2_L1.txt,nest3_L1.txt
     --ifo H1L1 --data nest1_H1L1.txt,nest2_H1L1.txt,nest3_H1L1.txt
  """
  parser.add_option("-d","--data",dest="data",action="append",help="A list " 
                    "of nested sampling data files for a particular IFO " 
                    "separated by commas (each IFO should use a separate " 
                    "instance of --data)")

  # number of live points
  parser.add_option("-N", "--Nlive", action="store", default=None, help="Number"
                    " of live points used in each parallel nested sampling "
                    "run. [required]", type="int", metavar="nlive",
                    dest="nlive")
  
  # Turn on 2D kdes
  #parser.add_option("--twodkdeplots", action="store_true", default=False,
  #                  dest="twodkdeplots")
  
  # Turn on R convergence tests
  #parser.add_option("--RconvergenceTests", action="store_true", default=False,
  #                  dest="RconvergenceTests")
  
  # get pulsar .par file
  parser.add_option("-p", "--parfile", dest="parfile", help="The TEMPO-style "
                    "pulsar parameter file used in the analysis. [required]", 
                    metavar="PNAME.par", default=None)
  
  # get pulsar correlation coefficient file
  parser.add_option("-c", "--corfile", dest="corfile", help="The TEMPO pulsar "
                    "correlation coefficient - if used in the analysis."
                    "[optional]", metavar="PNAME.mat", default=None)
  
  # get nested sampling analysis prior file
  parser.add_option("-P", "--priorfile", dest="priorfile", help="The prior "
                    "file used in the analysis. [required]",
                    metavar="prior.txt")
  
  # get heterodyned data files used for analysis
  parser.add_option("-B", "--Bkfiles", dest="Bkfiles", action="append", 
                    help="A heterodyned data file. [optional]", default=None)
  
  # get list of the detectors used (same order as Bk files)
  parser.add_option("-i", "--ifo", dest="ifos", help="The individual "
                    "interferometers from which analysis output, and data "
                    "files have been supplied. [required]", action="append",
                    default=None)
  
  # get number of bins for histogramming (default = 50)
  parser.add_option("-b", "--histbins", dest="histbins", help="The number of " \
                    "bins for histrogram plots. [default = 50]", default=50)
  
  # parse input options
  (opts, args) = parser.parse_args()
  
  # check that output path has been given
  if not opts.__dict__['outpath']:
    print "Must specify an output path"
    parser.print_help()
    exit(-1)
  
  # check that some data has been given
  if not opts.__dict__['data']:
    print "Must specify nested sampling data files"
    parser.print_help()
    exit(-1)
    
  # check that number of live points has been set
  if not opts.__dict__['nlive']:
    print "Must specify number of live points used in analysis"
    parser.print_help()
    exit(-1)
    
  # check that parfile has been given
  if not opts.__dict__['parfile']:
    print "Must specify a pulsar TEMPO .par file"
    parser.print_help()
    exit(-1)
    
  # check that a prior file has been given
  if not opts.__dict__['priorfile']:
    print "Must specify a prior file used during nested sampling"
    parser.print_help()
    exit(-1)
    
  if not opts.__dict__['ifos']:
    print "Must specify the interferometers analysed"
    parser.print_help()
    exit(-1)
  
  # upper limits wanted for output
  upperlimit=0.95
  
  # check that number of ifos is the same as the number of data lists
  nifos = len(opts.ifos)
  ndata = len(opts.data)
  
  if nifos != ndata:
    print "Number of IFOs and data lists are not equal"
    exit(-1)
  
  # sort out data into lists for each IFO
  data = []
  for i, ifo in enumerate(opts.ifos):
    data.append(opts.data[i].split(','))
  
  # call pulsarBayesPostProc function
  pulsarBayesPostProc( opts.outpath, data, \
                       upperlimit, opts.histbins, opts.priorfile, \
                       opts.parfile, opts.nlive, opts.Bkfiles, opts.ifos, \
                       opts.corfile ) #, \
                       #twodkdeplots=opts.twodkdeplots, \
                       #RconvergenceTests=opts.RconvergenceTests )
  