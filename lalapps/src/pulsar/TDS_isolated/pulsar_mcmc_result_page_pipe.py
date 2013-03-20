#!/usr/bin/env python

# -*- coding: utf-8 -*-
#
#       pulsar_mcmc_results_page_pipe.py
#
#       Copyright 2013
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

# This script will create a dag to produce a results page from the known pulsar
# analysis (using the MCMC results)

# import required modules
import sys
import os
import getopt
import re
import string
import tempfile
import exceptions
import random
import math
import ConfigParser
from optparse import OptionParser

import numpy as np

from lalapps import pulsarpputils as pppu
from lalapps import pulsar_pipeline_utils as ppu
from glue import pipeline
from pylal import git_version

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

# main function
if __name__=='__main__':
  description = \
"""This script is for creating a Condor DAG to produce the results page for multiple known pulsars."""

  usage = "Usage: %prog [options]"

  parser = OptionParser( usage = usage, description = description,
                         version = __version__ )
  
  parser.add_option("-i", "--configfile", dest="configfile", help="The configuration "
                    "(.ini) file for the analysis.")
  
  # parse input options
  (opts, args) = parser.parse_args()
  
  if not opts.configfile:
    print >> sys.stderr, "Error, not .ini file given"
    sys.exit(1)
    
  inifile = opts.configfile
  
  # create the config parser object and read in the ini file
  try:
    # the 'allow_no_value' arg only appeared in python 2.7, so check this
    if sys.version_info < (2, 7):
      cp = ConfigParser.ConfigParser()
      print >> sys.stderr, "Python version is < 2.7, so make sure configuration\
 file has un-needed values left out."
    else:
      cp = ConfigParser.ConfigParser(allow_no_value=True)
    cp.read(inifile)
  except:
    print >> sys.stderr, "Cannot read configuration file %s" % inifile
    sys.exit(1)
  
  """
  get the general inputs
  """
  try:
    # the output directory
    outdir = cp.get('general', 'outdir')
  except:
    outdir = None
  
  if not outdir:
    print >> sys.stderr, "No output directory specified!"
    sys.exit(1)
  
  try:
    # a directory of, or an individual, pulsar parameter file(s)
    parfile = cp.get('general', 'parfile')
  except:
    parfile = None
  
  if not parfile:
    print >> sys.stderr, "No parameter file/directory specified!"
    sys.exit(1)
  
  """
  log files
  """
  try:
    # a dag log directory basename
    logdir = cp.get('log', 'logdir')
  except:
    logdir = None
    
  if not logdir:
    print >> sys.stderr, "No DAG log directory specified!"
    sys.exit(1)
  
  """
  get individual pulsar results page-specific values 
  """
  try:
    # the results page creation code
    rexec = cp.get('resultspages', 'exec')
  except:
    rexec = None
    
  if not rexec or not os.path.isfile(rexec):
    print >> sys.stderr, "No results page executive script given!"
    sys.exit(1)
  
  try:
    # a list of MCMC directories (separated by commas)
    # mcmcdirs = /home/matthew/mcmc1, /home/matthew/mcmc2
    mcmcdirsall = cp.get('resultspages', 'mcmcdirs')
    mcmcdirs = (re.sub(' ', '', mcmcdirsall)).split(',')
  except:
    mcmcdirs = None
  
  if not mcmcdirs:
    print >> sys.stderr, "No MCMC file input directories specified!"
    sys.exit(1)
  
  try:
    # a list of interferometers (separated by commas)
    ifosall = cp.get('resultspages', 'ifos')
    ifosr = (re.sub(' ', '', ifosall)).split(',')
  except:
    ifosr = []
  
  if not ifosr:
    print >> sys.stderr, "No interferometers specified!"
    sys.exit(1)
  
  try:
    # a list of Bk directories for each IFO (in the same order as the IFOs are listed)
    bkdirsall = cp.get('resultspages', 'bkdirs')
    bkdirs = (re.sub(' ', '', bkdirsall)).split(',')
  except:
    bkdirs = None
  
  if len(bkdirs) != len(ifosr):
    print >> sys.stderr, "Heterodyned data directories and IFOs not consistent!"
    sys.exit(1)

  try:
    # a directory containing prior parameter files
    priordir = cp.get('resultspages', 'priordir')
  except:
    priordir = None

  try:
    # number of histogram bins
    histbins = cp.get('resultspages', 'histbins')
  except:
    histbins = None
    
  try:
    # using software injections
    swinj = cp.getboolean('resultspages', 'swinj')
  except:
    swinj = False
    
  try:
    # using hardware injections
    hwinj = cp.getboolean('resultspages', 'hwinj')
  except:
    hwinj = False
    
  try:
    # output figures in eps format as well as png
    epsoutr = cp.getboolean('resultspages', 'epsout')
  except:
    epsoutr = False
    
  """
  get collated results page-specific values
  """
  try:
    # the collate page creation code
    cexec = cp.get('collatepage', 'exec')
  except:
    cexec = None
    
  if not cexec or not os.path.isfile(cexec):
    print >> sys.stderr, "No collation page executive script given!"
    sys.exit(1)
    
  try:
    # get the results table order sorting format
    sorttype = cp.get('collatepage', 'sorttype')
  except:
    sorttype = None
    
  try:
    # get whether to compile the LaTeX results table
    compilelatex = cp.getboolean('collatepage', 'compilelatex')
  except:
    compilelatex = False
    
  try:
    # get the IFO(s) that you want to output a results table for (separated by commas)
    ifolistall = cp.get('collatepage', 'ifos')
    ifolist = (re.sub(' ', '', ifolistall)).split(',')
  except:
    ifolist = []
    
  try:
    # get a list of the upper limit values that you want output (comma separated)
    limlistall = cp.get('collatepage', 'outputlims')
    limlist = (re.sub(' ', '', limlistall)).split(',')
  except:
    limlist = []
    
  try:
    # get a list of the output pulsar values that you want output (comma separated)
    valslistall = cp.get('collatepage', 'outputvals')
    valslist = (re.sub(' ', '', valslistall)).split(',')
  except:
    valslist = []
    
  try:
    # output plots of the upper limits in histogram form
    histplot = cp.getboolean('collatepage', 'histplot')
  except:
    histplot = False
    
  try:
    # output plots of the h0 upper limits
    ulplot = cp.getboolean('collatepage', 'ulplot')
  except:
    ulplot = False
    
  try:
    # say whether to output prior limits on the above plots
    plotprior = cp.getboolean('collatepage', 'plotprior')
  except:
    plotprior = False
    
  try:
    # say whether to output plot in eps format as well as png
    epsout = cp.getboolean('collatepage', 'epsout')
  except:
    epsout = False
  
  # create the dag
  # create a log file that the Condor jobs will write to
  basename = re.sub(r'\.ini',r'',inifile) # remove .ini from config_file name

  if os.path.isdir(logdir):
    logfile = logdir + '/' + basename + '.dag.log'
    fh = open(logfile, "w" ) # creates file
    fh.close()
  else:
    print >> sys.stderr, "DAG log directory does not exist!"
    sys.exit(1)

  # create the DAG writing the log to the specified directory
  dag = pipeline.CondorDAG(logfile)
  dag.set_dag_file(basename)
  
  # set up jobs
  resultsjob = ppu.createresultspageJob(rexec, logdir)
  collatejob = ppu.collateresultsJob(cexec, logdir)
  
  inlimlist = False
  for lim in limlist:
    if lim in ['h95', 'ell' ,'sdrat', 'q22', 'sdpow']:
      inlimlist = True
      continue
    else:
      print >> sys.stderr, "Limit %s not recognised!" % lim
      inlimlist = False
      break
  
  if inlimlist:
    collatejob.add_arg('$(macrou)')
     
  invallist = False
  for val in valslist:
    if val in ['f0rot', 'f0gw', 'f1rot', 'f1gw', 'dist', 'sdlim', 'ra', 'dec', 'h0prior']:
      invallist = True
      continue
    else:
      print >> sys.stderr, "Parameter value %s not recognised!" % val
      invallist = False
      break
  
  if invallist:
    collatejob.add_arg('$(macron)')
  
  collatepage = ppu.collateresultsNode(collatejob)
  
  # check par file directory
  param_files = []
  if os.path.isfile(parfile):
    param_files.append(parfile)
  else:
    if os.path.isdir(parfile):
      # list the parameter files in the directory
      pfiles = os.listdir(parfile)

      for f in pfiles:
        if '.par' in f: # only use .par files
          param_files.append(os.path.join(parfile, f))
    else:
      print >> sys.stderr, "No par file or directory specified!"
      sys.exit(1)
  
  # set collate page values
  collatepage.set_parfile(parfile)
  
  if not os.path.isdir(outdir):
    print >> sys.stderr, "Output directory %s does not exist!" % outdir
    sys.exit(1)
  else:
    collatepage.set_outpath(outdir)
    collatepage.set_inpath(outdir)
    
  if sorttype:
    if sorttype in ['name', 'freq', 'ra', 'dec']:
      collatepage.set_sorttype(sorttype)
    else:
      print >> sys.stderr, "Invalid sort type %s - using default" % sorttype
    
  for ifo in ifolist:
    # check ifos are in the individual results page list 
    if ifo != 'Joint': # ifosr should not contain 'Joint'
      if ifo not in ifosr:
        print >> sys.stderr, "%d is not in the results page IFO list" % ifo
        sys.exit(1)

  collatepage.set_ifos(ifolist)
  
  if inlimlist:
    collatepage.set_outputlims(limlist)
  
  if invallist:
    collatepage.set_outputvals(valslist)

  if compilelatex:
    collatepage.set_compilelatex()
  
  if histplot:
    collatepage.set_outputhist()
    
  if ulplot:
    collatepage.set_outputulplot()
  
  if plotprior:
    collatepage.set_withprior()
    
  if epsout:
    collatepage.set_epsout()
  
  # loop over pulsars in par file directory
  for pfile in param_files:
    resultsnode = None
    resultsnode = ppu.createresultspageNode(resultsjob)
    
    resultsnode.set_parfile(pfile)
    resultsnode.set_outpath(outdir)
    
    for mcmcdir in mcmcdirs:
      if not os.path.isdir(mcmcdir):
        print >> sys.stderr, "MCMC directory %s not found!" % mcmcdir
        sys.exit(1)

    resultsnode.set_mcmcdir(mcmcdirs)
    
    for bkdir in bkdirs:
      if not os.path.isdir(bkdir):
        print >> sys.stderr, "Fine heterodyne directory %s not found!" % bkdir
        sys.exit(1)

    resultsnode.set_bkfiles(bkdirs)
    
    if priordir and os.path.isdir(priordir):      
      resultsnode.set_priordir(priordir) 
    
    resultsnode.set_ifos(ifosr)
    
    if histbins:
      resultsnode.set_histbins(histbins)
      
    if epsoutr:
      resultsnode.set_epsout()
      
    if hwinj:
      resultsnode.set_hwinj()
      
    if swinj:
      resultsnode.set_swinj()
    
    collatepage.add_parent(resultsnode)
    
    dag.add_node(resultsnode)
    
  dag.add_node(collatepage)
  
  # write out DAG
  dag.write_sub_files()
  dag.write_dag()
