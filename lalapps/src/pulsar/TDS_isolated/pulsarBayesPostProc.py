#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       cbcBayesPostProc.py
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
import numpy

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt

#local application/library specific imports
from pylal import SimInspiralUtils
from pylal import bayespputils as bppu
from pylal import git_version

import pulsarpputils as pppu

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date


def pulsarBayesPostProc( outdir, data, oneDMenu, twoDGreedyMenu, GreedyRes,
                         confidence_levels, twoDplots,
                         #nested sampling options
                         ns_Nlive=None,
                         #Turn on 2D kdes
                         twodkdeplots=False,
                         #Turn on R convergence tests
                         RconvergenceTests=False,
                         # analysis data
                         Bkdata=None, ifos=None
                         # pulsar information
                         parfile=None, corfile=None
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
  
  # if input format is xml
  votfile = None
  if '.xml' in data[0]:
    peparser = bppu.PEOutputParser('xml')
    commonResultsObj = peparser.parse(data[0])
    thefile = open(data[0],'r')
    votfile = thefile.read()
  else: # input as standard nested sampling input
    peparser = bppu.PEOutputParser('ns')
    commonResultsObj = peparser.parse(data, Nlive=ns_Nlive, xflag=False)

  # Create an instance of the posterior class using the posterior values loaded
  # from the file.
  pos = bppu.Posterior( commonResultsObj, SimInspiralTableEntry=None, \
                        votfile=votfile )

  # convert phi0' and psi' to phi0 and psi
  
  # read in par file if given
  if parfile is not None:
    par = pppu.psr_par(parfile)
    
    # get par file info as a string as well (might need!)
    parstr = str(par)
  
  # convert phi0' and psi' to phi0 and psi

    
    
# main function
if __name__=='__main__':
  from optparse import OptionParser
  
  parser=OptionParser()
  parser.add_option("-o", "--outpath", dest="outpath", help="Make page and " \
                    "plots in DIR", metavar="DIR")
  parser.add_option("-d","--data",dest="data",action="append",help="A list " \
                    "of nested sampling data files")
  
  parser.add_option("--no2D", action="store_true", default=False, help="Skip " \
                    "2-D plotting.")

  # number of live points
  parser.add_option("--Nlive", action="store", default=None, help="Number of " \
                    "live points used in each parallel nested sampling run.", \
                    type="int")
  
  # Turn on 2D kdes
  parser.add_option("--twodkdeplots", action="store_true", default=False,
                    dest="twodkdeplots")
  
  # Turn on R convergence tests
  parser.add_option("--RconvergenceTests", action="store_true", default=False,
                    dest="RconvergenceTests")
  
  # get pulsar .par file
  parser.add_option("-p", "--parfile", dest="parfile", help="The pulsar " \
                    "parameter file.", default=None)
  
  # get pulsar correlation coefficient file
  parser.add_option("-c", "--corfile", dest="corfile", help="The pulsar " \
                    "correlation coefficient file.", default=None)
  
  # get heterodyned data files used for analysis
  parser.add_option("-B", "--Bkfiles", dest="Bkfiles", action="append", 
                    help="A heterodyned data file.", default=None)
  
  # get list of the detectors used (same order as Bk files)
  parser.add_option("-i", "--ifos", dest="ifos", action="append", default=None)
  
  # parse input options
  (opts,args)=parser.parse_args()
  
  # Confidence levels wanted for output
  confidenceLevels=[0.67,0.9,0.95,0.99]
  
  # call pulsarBayesPostProc function
  cbcBayesPostProc( opts.outpath, opts.data, oneDMenu, twoDGreedyMenu,
                    greedyBinSizes, confidenceLevels, twoDplots,
                    ns_Nlive=opts.Nlive,
                    twodkdeplots=opts.twodkdeplots,
                    RconvergenceTests=opts.RconvergenceTests,
                    opts.Bkfiles, opts.ifos,
                    opts.parfile, opts.corfile )
  