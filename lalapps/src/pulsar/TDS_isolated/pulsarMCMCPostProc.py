#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       pulsarMCMCPostProc.py
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
import re
import shutil
import math
import datetime

#related third party imports
import numpy as np
import subprocess as sp

#local application/library specific imports
from pylal import bayespputils as bppu
from pylal import git_version

from lalapps import pulsarpputils as pppu

import urllib2
from BeautifulSoup import BeautifulSoup as bs

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

# create css file text
csstext = """
/* create body style */
body {
  font-family: Verdana, Geneva, "Trebuchet MS", sans-serif;
}

/* create header name style */
h1 {
  margin: 0 0 16px 0;
  padding: 0 0 16px 0;
  font-size: 20px;
  font-weight: bold;
  letter-spacing: 0px;
  border-bottom: 1px solid #999;
  font-family: Verdana, Geneva, "Trebuchet MS", Sans-Serif;
  text-shadow: 2px 2px 2px #ccc;
}

h2 {
  margin: 0 0 12px 0;
  padding: 0 0 12px 0;
  font-size: 16px;
  font-weight: bold;
  border-bottom: 1px solid #999;
  font-family: Verdana, Geneva, "Trebuchet MS", Sans-Serif;
  text-shadow: 1px 1px 1px #ccc;
}

# create pulsar parameter display style
/* .par {
  float: left;
  padding: 10px;
} */

/* .limits {
  float: left;
} */

.wrapper {
 /* float: left; */
}

/* create footer style */
#footer {
  border-top: 1px solid #999;
  padding: 15px;
  font-family: monospace;
  text-align: center;
}
  
/* create link style */
a:link{
  color: #000000;
  text-decoration: none;
}

a:visited{
  color: #000000;
  text-decoration: none;
}

a:hover{
  color: #000000;
  text-decoration: none;
  text-shadow: 2px 2px 2px #ccc;
}

/* create a class for a posterior plot image */
.posplot{
  height: 275px;
  /* padding: 8px; */
  border: 0px solid #999;
  /* box-shadow: 2px 2px 2px #888; */
}
  
/* create a class for a rotated image */
.rotated90clockwise{
  padding: 0px;
  border: 0px solid #999;
  /* box-shadow: 2px 2px 2px #888; */
  width: 275px;
  -moz-transform: rotate(90deg);
  -webkit-transform: rotate(90deg);
  -o-transform: rotate(90deg);
  -ms-transform: rotate(90deg);
}

.rotated90anticlockwise{
  padding: 0px;
  width: 275px;
  -moz-transform: rotate(270deg);
  -webkit-transform: rotate(270deg);
  -o-transform: rotate(270deg);
  -ms-transform: rotate(270deg);
}

/* create a class for a Bk data plot */
.dataplot{
  height: 275px;
}

/* create class for an amplitude spectral density plot */
.asdplot{
  height: 275px;
}

/* create class for an MCMC chain plot */
.chainplot{
  height: 275px;
}

/* set defaults for table data and table headers */
table{
  border: 0px;
}

td{
  padding: 0px 8px 0px 8px;
  border: 0px;
}

th{
  padding: 0px 8px 0px 8px;
  border: 0px;
}

.leftborder{
  border-left: 1px solid #000;
}

.rightborder{
  border-right: 1px solid #000;
}

.topborder{
  border-top: 1px solid #000;
}

.bottomborder{
  border-bottom: 1px solid #000;
}
  
/* set text colour classes for detectors */
.H1{
  color: red;
}

.H2{
  color: cyan;
}

.L1{
  color: green;
}

.V1{
  color: blue;
}

.G1{
  color: magenta;
}

.Joint{
  color: black;
  font-weight: bold;
}
"""

# convert a floating point number into a string in X.X x 10^Z format
def exp_str(f, p=1):
  if p > 16:
    print >> sys.stderr, "Precision must be less than 16 d.p."
    p = 16
  
  s = '%.16e' % f
  ssplit = s.split('e')
  return '%.*f&times;10<sup>%d</sup>' % (p, float(ssplit[0]), int(ssplit[1]))

def exp_latex_str(f, p=1):
  if p > 16:
    print >> sys.stderr, "Precision must be less than 16 d.p."
    p = 16
  
  s = '%.16e' % f
  ssplit = s.split('e')
  return '\\ensuremath{%.*f\!\\times\!10^{%d}}'  % (p, float(ssplit[0]), \
int(ssplit[1]))
  
# convert a right ascension string in format 'hh:mm:ss.s' to a html string like
# H^h M^m S^s.ss
def ra_htmlstr(ra):
  hms = ra.split(":")
  
  if len(hms) == 1:
    hms.append('0')
    hms.append('0')
  elif len(hms) == 2:
    hms.append('0')

  ss = ('%.2f' % float(hms[2])).split('.')
  
  return "%s<sup>h</sup>%s<sup>m</sup>%s<sup>s</sup>.%s" % (hms[0].zfill(2), \
hms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2) )

# convert a declination string in format 'dd:mm:ss.s' to a html string like
# dd^o mm' ss''.ss
def dec_htmlstr(ra):
  dms = ra.split(":")
  
  if len(dms) == 1:
    dms.append('0')
    dms.append('0')
  elif len(dms) == 2:
    dms.append('0')

  ss = ('%.2f' % float(dms[2])).split('.')
    
  return "%s&deg;%s'%s\".%s" % ((re.sub('\+', '', dms[0])).zfill(2), \
dms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2))


# return the correlation coefficient of all parameters (but make sure to not
# include the logl or post values. It takes in a posterior object.
def corrcoef(pos):
  header_string=''
  posterior_table=[]
 
  for param_name, one_pos in pos:
    if param_name != 'logl' and param_name != 'post':
      column = np.array(one_pos.samples)
      header_string += param_name + ' '
      posterior_table.append(column)
  
  posterior_table = tuple(posterior_table)
  parray = np.column_stack(posterior_table)
  
  cc = np.corrcoef(parray, rowvar=0, bias=1)
  
  return cc, header_string

  
# a function that attempt to load an MCMC chain file: first it tries using 
# numpy.loadtxt; if it fails it tries reading line by line and checking
# for consistent line numbers, skipping lines that are inconsistent; if this
# fails it returns None
def readmcmcfile(cf):
  cfdata = None
  
  # first try reading in with loadtxt (skipping lines starting with %)
  try:
    cfdata = np.loadtxt(cf, comments='%')
  except:
    try:
      fc = open(cf, 'r')
      
      # read in header lines and count how many values each line should have
      headers = fc.readline()
      headers = fc.readline() # column names are on the second line
      # remove % from start
      headers = re.sub('%', '', headers)
      # remove rads
      headers = re.sub('rads', '', headers)
      # remove other brackets e.g. around (iota)
      headers = re.sub('[()]', '', headers)
      
      lh = len(headers.split())
      
      cfdata = np.array([])
      
      lines = cf.readlines()
      
      for i, line in enumerate(lines):
        if '%' in line: # skip lines containing %
          continue
        
        lvals = line.split()
        
        # skip line if number of values isn't consistent with header
        if len(lvals) != lh:
          continue
        
        # convert values to floats
        try:
          lvalsf = map(float, lvals)
        except:
          continue
        
        # add values to array
        if i==0:
          cfdata = np.array(lvalsf)
        else:
          cfdata = np.vstack((cfdata, lvalsf))
      
      if cfdata.size == 0:
        cfdata = None
    except:
      cfdata = None

  return cfdata
          
      
# list of parameters to display (in this order)
paramdisplist = ['RAJ', 'DECJ', 'F0', 'F1', 'F2', 'PEPOCH', 'X' 'E' \
'EPS1', 'EPS2', 'OM', 'T0', 'TASC', 'PB']

# html text to display for different parameter names
paramtextdisp = {'RAJ': '&alpha;', 'DECJ': '&delta;', \
                 'RA': '&alpha;', 'DEC': '&delta;', \
                 'F0': 'f<sub>0</sub> (Hz)', 'F1': 'f<sub>1</sub> (H/s)', \
                 'F2': 'f<sub>2</sub> (H/s<sup>2</sup>)', \
                 'PEPOCH': 'epoch (MJD)', 'A1': 'a sin<it>i</i> (lt s)', \
                 'E': '<it>e</it>', 'EPS1': '&epsilon;<sub>1</sub>', \
                 'EPS2': '&epsilon;<sub>2</sub>', \
                 'T0': 'T<sub>0</sub> (MJD)', \
                 'TASC': 'T<sub>asc</sub> (MJD)', \
                 'OM': '&omega;<sub>0</sub>$deg;',
                 'PB': 'Period (days)', 'H0': 'h<sub>0</sub>', \
                 'COSIOTA': 'cos&iota;', 'PSI': '&psi; (rad)', \
                 'PHI0': '&phi;<sub>0</sub> (rad)', \
                 'PMRA': 'p.m. &alpha; (rad/s)', \
                 'PMDC': 'p.m. &delta; (rad/s)', \
                 'PMDEC': 'p.m. &delta; (rad/s)'}

# a class containing function to output parameter vales in the appropriate
# format
class paramdisp:
  def RAJ(f): return ra_htmlstr(f)
  def RA(f): return ra_htmlstr(f)
  def DECJ(f): return dec_htmlstr(f)
  def DEC(f): return dec_htmlstr(f)
  def F0(f): return '%.5f' % float(f)
  def F1(f): return exp_str(float(f), 2)
  def F2(f): return exp_str(float(f), 2)
  def PEPOCH(f): return '%d' % int(float(f)) # return epoch as an integer
  def A1(f): return '%.2f' % float(f)
  def E(f): return '%.2f' % float(f)
  def EPS1(f): return exp_str(float(f), 2)
  def EPS2(f): return exp_str(float(f), 2)
  def T0(f): return '%.3f' % float(f)
  def TASC(f): return '%.3f' % float(f)
  def OM(f): return '%.1f' % float(f)
  def PB(f): return '%.3f' % float(f)
  def H0(f): return exp_str(float(f), 2)
  def COSIOTA(f): return '%.2f' % float(f)
  def PHI0(f): return '%.2f' % float(f)
  def PSI(f): return '%.2f' % float(f)
  
# a class containing function to output parameter vales in the appropriate
# format (for use when outputting posterior stats)
class paramdisp2:
  def RAJ(f): return exp_str(float(f), 2)
  def DECJ(f): return exp_str(float(f), 2)
  def RA(f): return exp_str(float(f), 2)
  def DEC(f): return exp_str(float(f), 2)
  def F0(f): return exp_str(float(f), 2)
  def F1(f): return exp_str(float(f), 2)
  def F2(f): return exp_str(float(f), 2)
  def PEPOCH(f): return '%d' % int(float(f)) # return epoch as an integer
  def A1(f): return '%.2f' % float(f)
  def E(f): return '%.2f' % float(f)
  def EPS1(f): return exp_str(float(f), 2)
  def EPS2(f): return exp_str(float(f), 2)
  def T0(f): return '%.2f' % float(f)
  def TASC(f): return '%.2f' % float(f)
  def OM(f): return '%.2f' % float(f)
  def PB(f): return '%.2f' % float(f)
  def H0(f): return exp_str(float(f), 2)
  def COSIOTA(f): return '%.2f' % float(f)
  def PHI0(f): return '%.2f' % float(f)
  def PSI(f): return '%.2f' % float(f)
  def PMRA(f): return exp_str(float(f), 2)
  def PMDC(f): return exp_str(float(f), 2)
  def DEFAULT(f): return '%.2f' % float(f)
  
# concatenate a list of strings (default to a new line between each)
def cattext(strlist, delim='\n'):
  return delim.join([strpart for strpart in strlist]) 
  
# output figures - take in the figure, the output path, a filename string and a
# list out output types e.g. ['png', 'eps']. Return a dictionary of the
# filenames for the different types
def output_fig(myfig, outpath, fname, ftypes):
  fnameret = {}
  
  if not ftypes or not fname or not outpath:
    print >> sys.stderr, "Error, problem outputting figure"
    sys.exit(1)
  
  for i, t in enumerate(ftypes):
    fnameret[t] = fname + '.' + t
    plotpath = os.path.join(outpath, fnameret[t])
    
    try:
      myfig.savefig(plotpath)
    except:
      print >> sys.stderr, "Error outputting figure %s" % plotpath
      sys.exit(1)
  
  return fnameret
  
# main function
if __name__=='__main__':
  from optparse import OptionParser
  
  description = \
"""This script is for creating a results output page and plots for the known
   pulsar analysis. It uses inputs from the old MCMC code
   lalapps_pulsar_parameter_estimation."""

  epilog = \
"""
An example of usage for a case when three nested sampling runs have been
performed for two interferometers (H1 and L1):
  %s --ifo H1 --Bkfiles /home/me/hetdata/H1 --ifo L1 --Bkfiles
/home/me/hetdata/L1 --parfile /home/me/pardir --priorfile priors.txt --histbins
50 --mcmcdirs /home/me/MCMCchains/pdfs1 --mcmcdirs /home/me/MCMCchains/pdfs2
--outpath /home/me/public_html/results
""" % os.path.basename(sys.argv[0])
  usage = "Usage: %prog [options]"
  
  parser = OptionParser( usage = usage, description = description,
                         version = __version__, epilog = epilog )
  
  parser.add_option("-o", "--outpath", dest="outpath", help="The path for "
                    "the analysis output", metavar="DIR")
  
  """
   data will be read in in the following format:
     multiple MCMC directories from each IFO would be entered as follows:
     --ifo H1 --Bkfiles /home/user/finehetsH1
     --ifo L1 --Bkfiles /home/user/finehetsL1
     --mcmcdirs /home/user/pdfs1
     --mcmcdirs /home/user/pdfs2
     --mcmcdirs /home/user/pdfs3
  """
  parser.add_option("-m","--mcmcdirs", dest="mcmcdirs",\
                    action="append", help="A list "\
                    "of MCMC directories containing chain for all IFOs")
  
  # get pulsar .par file
  parser.add_option("-p", "--parfile", dest="parfile", help="The "
                    "directory containing several, or an individual, "
                    "TEMPO-style pulsar parameter file used in the analysis. "
                    "[required]", metavar="PNAME.par", default=None)
 
  # get heterodyned data files used for analysis
  parser.add_option("-B", "--Bkfiles", dest="Bkfiles", action="append", 
                    help="A directory containing the heterodyned data files "
                    "for a given detector [required]", default=None)
  
  # get list of the detectors used (same order as Bk files)
  parser.add_option("-i", "--ifo", dest="ifos", help="The individual "
                    "interferometers from which analysis output, and data "
                    "files have been supplied. [required]", action="append",
                    default=None)
  
  # get number of bins for histogramming (default = 50)
  parser.add_option("-b", "--histbins", dest="histbins", help="The number of " \
                    "bins for histrogram plots. [default = 50]", default=50)
  
  # say if you want to output eps figure (default to False)
  parser.add_option("-e", "--epsout", dest="epsout", help="Set if wanting to " \
                    "output eps versions of figures", action="store_true",
                    default=False)
  
  # say if the pulsars are software injections (default to False)
  parser.add_option("-s", "--sw-inj", dest="swinj", help="Set if the " \
                    "pulsars are software injections", action="store_true", \
                    default=False)
                    
  # say if the pulsars are hardware injections (default to False)
  parser.add_option("-w", "--hw-inj", dest="hwinj", help="Set if the " \
                    "pulsars are hardware injections", action="store_true", \
                    default=False)
  
  p_args=''
  for arg in sys.argv:
    p_args += arg + ' '
  
  # parse input options
  (opts, args) = parser.parse_args()
  
  # check that output path has been given
  if not opts.__dict__['outpath']:
    print "Must specify an output path"
    parser.print_help()
    exit(-1)
  else:
    outpath = opts.outpath
  
  # check that some data has been given
  if not opts.__dict__['mcmcdirs']:
    print "Must specify MCMC chain directories"
    parser.print_help()
    exit(-1)
    
  # check that parfile has been given
  if not opts.__dict__['parfile']:
    print "Must specify a pulsar TEMPO .par file"
    parser.print_help()
    exit(-1)
       
  if not opts.__dict__['ifos']:
    print "Must specify the interferometers analysed"
    parser.print_help()
    exit(-1)
  else:
    ifos = opts.ifos
  
  if not opts.__dict__['Bkfiles']:
    print "Must specify the heterodyned data files"
    parser.print_help()
    exit(-1)
  else:
    Bkfiles = opts.Bkfiles

  if opts.epsout:
    ftypes = ['png', 'eps']
  else:
    ftypes = ['png']
  
  swinj = opts.swinj
  hwinj = opts.hwinj
  
  # create output directory if it doesn't exist
  if not os.path.isdir(outpath):
    try:
      os.mkdir(outpath)
    except:
      print >> sys.stderr, "Cannot create output directory %s" % outpath
      sys.exit(1)
    
  # check that number of ifos is the same as the number of data lists
  nifos = len(ifos)
  ndata = len(Bkfiles)
  
  if nifos != ndata:
    print "Number of IFOs and data lists are not equal"
    exit(-1)
  
  ifosNew = list(ifos)
  if nifos > 1:
    ifosNew.append('Joint')
  
  # split MCMC directories
  mcmcdirs = opts.mcmcdirs
  
  # check if parfile is a single file or a directory
  param_files = []
  if os.path.isfile(opts.parfile):
    param_files.append(opts.parfile)
  else:
    if os.path.isdir(opts.parfile):
      # list the parameter files in the directory
      pfiles = os.listdir(opts.parfile)
      
      for f in pfiles:
        param_files.append(opts.parfile + '/' + f)
    else:
      print >> sys.stderr, "No par file or directory specified!"
      sys.exit(1)
  
  # create CSS
  cssname = 'psr.css'
  try:
    cssfile = os.path.join(outpath, cssname)
    css = open(cssfile, 'w')
  except:
    print >> sys.stderr, "Cannot open CSS file %s" % cssfile
    sys.exit(1)

  css.write(csstext);
  css.close()
  
  # get time/date for file creation
  now = datetime.datetime.now()
  
  # create html page and table for each result
  # create html page
  htmlpage = os.path.join(outpath, 'index.html')
  try:
    htmlout = open(htmlpage, 'w')
  except:
    print >> sys.stderr, 'Could not open html page %s' % htmlpage
    
  htmltext = []
  htmltext.append( \
"""
<!DOCTYPE html>
<html lang="en">

<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
  <meta name="description" content="Results page"/>
  <meta charset="UTF-8"> 

  <title>Results page</title>

  <link rel="stylesheet" type="text/css" href="%s"/>
</head>
""" % cssname )
    
  htmltext.append('<body>\n')
  htmltext.append('<table>\n')
  
  # create LaTeX table page of results
  latexpage = os.path.join(outpath, 'resultstable.tex')
  try:
    latexout = open(latexpage, 'w')
  except:
    print >> sys.stderr, 'Could not open LaTeX page %s' % htmlpage
  
  latextext = []
  latextext.append( \
"""
%\documentclass[preprint,12pt]{aastex}
\documentclass{emulateapj}

\usepackage{amsmath}
% The rotate command for deluxetable does not work in emulateapj, so use a 
% landscape mode instead
\usepackage{lscape}

\\begin{document}
""")
  
  # list of amplitude spectral densities
  asds = []
  ullist = []
  sdlist = []
  
  firstpsr = True
  
  # loop through par files to produce plots
  for parfile in param_files:
    # read data from par file
    try:
      par = pppu.psr_par(parfile)
    except:
      print >> sys.stderr, "Par file %s could not be opened!" % parfile
      continue # move on to next pulsar
   
    pname = par['PSRJ']
    if not pname:
      print >> sys.stderr, "No PSRJ value in par file %s" % parfile
      continue # move on to next pulsar
    
    # check if heterodyned data exists for this pulsar and each detector
    Bkdata = []
    nofile = False
    for i, ifo in enumerate(ifos):
      Bkdata.append(Bkfiles[i] + '/finehet_' + pname + '_' + ifo)
      
      # check files exist if not then skip the pulsar
      if not os.path.isfile(Bkdata[i]):
        nofile = True
        break
    if nofile:
      continue
    
    # check that MCMC chains exist for this pulsar and each detector (including
    # joint)
    for ifo in ifosNew:
      for chaindir in mcmcdirs:
        cfile = chaindir + '/' + 'MCMCchain_' + pname + '_' + ifo
        
        if not os.path.isfile(cfile):
          nofile = True
          break
    if nofile:
      continue
    
    # check required parameters
    f0 = par['F0']
    if not f0:
      print >> sys.stderr, "No F0 value in par file %s" % parfile
      continue # move on to next pulsar
      
    f1 = par['F1']
    if not f1:
      print >> sys.stderr, "No F1 value in par file %s" % parfile
      continue # move on to next pulsar
    
    # create output directory for pulsar
    puldir = os.path.join(outpath, pname)
    if not os.path.isdir(puldir):
      os.mkdir(puldir)
    
    # attempt to get pulsar distance and proper motion corrected age from the
    # ATNF catalogue
    ages = None
    dists = None
    try:
      atnfurl = \
'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?Dist=Dist&Age_i=\
Age_i&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=&\
sort_order=asc&condition=&pulsar_names=' + re.sub('\+', '%2B', pname) + \
'&ephemeris=selected&\
submit_ephemeris=Get+Ephemeris&coords_unit=raj%2Fdecj&radius=&coords_1=\
&coords_2=&style=Long+with+last+digit+error&no_value=*&fsize=3&x_axis=&x_scale=\
linear&y_axis=&y_scale=linear&state=query' 
      soup = bs(urllib2.urlopen(atnfurl).read())
      pdat = soup.pre # data exists in the pre html environment
      
      for line in pdat:
        vals = line.split('\n') # split at any new lines
        
        for row in vals:
          if 'DIST' in row:
            dists = row.split() # split row at whitespace
            
          if 'AGE_I' in row:
            ages = row.split()
    except:
      print >> sys.stderr, 'Problem accessing ATNF!'
      sys.exit(1)
    
    if ages:
      if len(ages) > 1:
        try:
          age = float(ages[1])
        except:
          print >> sys.stderr, "Age not a number!"
          sys.exit(1)
        f1sd = -f0/(2*age*365.25*86400)
      else:
        f1sd = f1
    # check if f1 is positive - in which case set based on conservative age of
    # tau = f0/(2*f1) = 10^9 years (for GC pulsars)
    elif f1 > 0 and pname != 'J1824-2452A':
      f1sd = -f0/(2*1e9*365.25*86400) # set f1 for spin-down limi calc
    else:
      f1sd = f1
   
    # first try getting a distance from a par file
    dist = par['DIST']
    # if not available look in download value
    if not dist:
      # see if it was in the ATNF catalogue
      if dists:
        if len(dists) > 1:
          try:
            dist = float(dists[1])
          except:
            print >> sys.stderr, "Distance not a number!"
            sys.exit(1)
        else:
          print >> sys.stderr, "Distance not a number!"
          sys.exit(1)
    
    # if on first available pulsars create the output page table header
    if firstpsr:
      htmltext.append('<tr>\n<th class="bottomborder" colspan="5"></th>')
      latextext.append( \
"""\
\clearpage
\LongTables %% table can span multiple pages
\\begin{landscape}
\\begin{deluxetable}{%s}
%%\\rotate
\\tabletypesize{\\footnotesize}
\\tablewidth{0pt}
\\tablecaption{Limits on the gravitational wave amplitude for known pulsars}
\\tablehead{\multicolumn{5}{c}{~} \
""" % ' '.join(['c']*(5+3*len(ifosNew))) )

      for ifo in ifosNew:
        htmltext.append( \
"""<th class="%s" colspan="3" style="text-align:center; border-left:1px \
solid #000; border-bottom:1px solid #000">%s</th>""" % (ifo, ifo))
        latextext.append(' & \multicolumn{3}{c}{%s}' % ifo)

      htmltext.append( \
"""
  </tr>
  <tr>
    <th class="bottomborder">Pulsar</th>
    <th class="bottomborder">Frequency (Hz)</th>
    <th class="bottomborder">Spin-down (Hz/s)</th>
    <th class="bottomborder">Distance (kpc)</th>
    <th class="bottomborder">Spin-down limit</th>
""" )
      latextext.append( \
"""\\\\
\colhead{Pulsar} & \colhead{$\\nu$ (Hz)} & \colhead{$\dot{\\nu}$ (Hz/s)} & \
\colhead{distance (kpc)} & \colhead{spin-down limit} \
""" )
      
      for ifo in ifosNew:
        htmltext.append( \
"""
  <th class="leftborder" style="border-bottom:1px solid \
#000">h<sub>0</sub><sup>95%</sup></th>
  <th class="bottomborder">&#949;</th>
  <th class="bottomborder">ratio</th>
""" )
        latextext.append( \
"""\
& \colhead{$h_0^{95\%}$} & \colhead{$\\varepsilon$} & \colhead{ratio} \
""" )

      htmltext.append('</tr>')
      latextext.append('}\n\startdata')
      
      firstpsr = False
    
    # create html page
    htmlppage = os.path.join(puldir, 'index.html')
    try:
      htmlpout = open(htmlppage, 'w')
    except:
      print >> sys.stderr, 'Could not open html page %s' % htmlppage
    
    htmlptext = []
    
    htmlptext.append( \
"""
<!DOCTYPE html>
<html lang="en">

<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
  <meta name="description" content="PSR %s"/>
  <meta charset="UTF-8"> 

  <title>PSR %s</title>

  <link rel="stylesheet" type="text/css" href="../%s"/>
 
</head>""" % (pname, pname, cssname))
    
    htmlptext.append('<body>')
    
    # add in javascript to allow text toggling
    htmlptext.append( \
"""
<script language="javascript"> 
function toggle(id) {
  var ele = document.getElementById(id);
  if(ele.style.display == "block"){
    ele.style.display = "none";
  }
  else
    ele.style.display = "block";
}
</script>
""" )
    
    # print out pulsar name to file
    pheadertext = []
    if not swinj and not hwinj:
      pulsaratnflink = \
'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?startUserDefined=\
true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&\
pulsar_names=' + re.sub('\+', '%2B', pname) + \
'&ephemeris=long&submit_ephemeris=Get+Ephemeris&\
coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2=&style=Long+with+last+digit+\
error&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=\
query'
    
    # set title prefix depending is injection or not
    psrnameprefix = 'PSR'
    if swinj:
      psrnameprefix = 'SWINJ'
    elif hwinj:
      psrnameprefix = 'HWINJ'

    if not swinj and not hwinj:
      pheadertext.append('<h1><a href="%s">%s %s</a></h1>\n' %
(pulsaratnflink, psrnameprefix, pname))
    else:
      pheadertext.append('<h1>%s %s</h1>\n' % (psrnameprefix, pname))
    
    # copy par file to pulsar directory
    shutil.copy(parfile, os.path.join(puldir, pname + '.par'))

    # print out par file info
    partext = []
    partext.append('<a href="'+pname+'.par">')
    partext.append('<table>')
    
    for param in paramdisplist:
      pa = par[param]
      if pa:
        partext.append('<tr>\n<td>%s</td>' % paramtextdisp[param])
        dispfunc = paramdisp.__dict__[param]
        partext.append('<td>%s</td>\n</tr>' % dispfunc(str(pa)))
    partext.append('</table></a>')
    
    # get spin-down upper limit
    sdlim = pppu.spin_down_limit(f0, f1sd, dist)
    sdlist.append(sdlim)
    
    # print out info to the main table
    iscor = ''
    iscorlatex = ''
    if f1sd != f1:
      iscor = '<sup>&dagger;</sup>'
      iscorlatex = '^{\dagger}'
    
    htmltext.append( \
"""
<tr>
  <td><a href="%s/">%s</a></td>
  <td>%.2f</td>
  <td>%s</td>
  <td>%.1f</td>
  <td>%s%s</td>
""" % (pname, pname, f0, exp_str(f1,2), dist, exp_str(sdlim), iscor) )

    latextext.append( \
"""\
%s & $%.2f$ & $%s$ & $%.1f$ & $%s$$%s$ \
""" % (re.sub('\-', '\\textminus', pname), f0, exp_latex_str(f1, 2), dist, \
exp_latex_str(sdlim), iscorlatex) )
    
    # get time series and PSD plots
    Bkdata = []
    plotpsds = True
    plotfscan = True
    
    for i, ifo in enumerate(ifos):
      Bkdata.append(Bkfiles[i] + '/finehet_' + pname + '_' + ifo)
    asdtime = 14400 # set time over which to produce the asds
    
    Bkfigs, psdfigs, fscanfigs, asdlist = pppu.plot_Bks_ASDs( Bkdata, ifos, \
asdtime, plotpsds=plotpsds, plotfscan=plotfscan, removeoutlier=8 )
    
    if asdlist:
      asds.append(asdlist)
    
    # output plots of time series and psd
    Bkfigname = []
    psdfigname = []
    fscanfigname = []
    for i, ifo in enumerate(ifos):
      figname = output_fig(Bkfigs[i], puldir, 'Bk_'+ifo, ftypes)
      Bkfigname.append(figname)
      
      if plotpsds:
        figname = output_fig(psdfigs[i], puldir, 'ASD_'+ifo, ftypes)
        psdfigname.append(figname)
        
        if plotfscan:
          figname = output_fig(fscanfigs[i], puldir, 'fscan_'+ifo, ftypes)
          fscanfigname.append(figname)
    
    # loop over detectors
    poslist = []
    mcmcgr = [] # list of Gelman-Rubins stats for detector
    nefflist = [] # effective sample size for each combined chain
    chainlens = [] # original average length of each chain
    
    # get the posterior statistics i.e. maxL, means, std. devs.
    max_pos = []
    maxLvals = []
    meanvals = []
    stddevvals = []
    medvals = []
    corrcoefs = []
    corrcoefsheader = []
    
    erridx = False
    
    for i, ifo in enumerate(ifosNew):
      # read in MCMC chains
      mcmc = []
      neffs = [] # number of effective samples for each chain
      cl = []
      grr = {}
      for chaindir in mcmcdirs:
        cfile = chaindir + '/' + 'MCMCchain_' + pname + '_' + ifo
        
        if os.path.isfile(cfile):
          # load MCMC chain
          mcmcChain = readmcmcfile(cfile)
          
          # if no chain was returned write out an error message and skip this
          # pulsar
          if mcmcChain == None:
            outerrfile = os.path.join(outpath, 'mcmc_chain_errors.txt')
            errfile = open(outerrfile, 'a') # append to file
            errfile.write('MCMC file %s could not be read in\n' % cfile)
            errfile.close()
            erridx = True
            break
          
          # find number of effective samples for the chain
          neffstmp = []
          for j in range(1, mcmcChain.shape[1]):
            neff, acl, acf = bppu.effectiveSampleSize(mcmcChain[:,j])
            neffstmp.append(neff)
            
          # get the minimum effective sample size
          neffs.append(min(neffstmp))
          
          nskip = math.ceil(mcmcChain.shape[0]/min(neffstmp))
          
          # output every nskip (independent) value
          mcmc.append(mcmcChain[::nskip,:])
          cl.append(mcmcChain.shape[0])
        else:
          print >> sys.stderr, "File %s does not exist!" % cfile
          sys.exit(1)
      
      if erridx:
        break
      
      nefflist.append(math.fsum(neffs))
      chainlens.append(np.mean(cl))
      
      # output data to common results format
      # get first line of MCMC chain file for header names
      cf = open(cfile, 'r')
      headers = cf.readline()
      headers = cf.readline() # column names are on the second line
      # remove % from start
      headers = re.sub('%', '', headers)
      # remove rads
      headers = re.sub('rads', '', headers)
      # remove other brackets e.g. around (iota)
      headers = re.sub('[()]', '', headers)
      cf.close()
      
      # get Gelman-Rubins stat for each parameter
      for idx, parv in enumerate(headers.split()):
        lgr = []
        if parv != 'logL':          
          for j in range(0, len(mcmc)):
            achain = mcmc[j]
            singlechain = achain[:,idx]
            lgr.append(singlechain)
          grr[parv.lower()] = pppu.gelman_rubins(lgr)
      
      mcmcgr.append(grr)
      
      # logL in chain is actually log posterior, so also output the posterior
      # values (can be used to estimate the evidence)
      headers = headers.replace('\n', '\tpost\n')
      
      # output full data to common format
      comfile = os.path.join(puldir, 'common_tmp.dat')
      try:
        cf = open(comfile, 'w')
      except:
        print >> sys.stderr, "Can't open commom posterior file!"
        sys.exit(1)
          
      cf.write(headers)
      for narr in mcmc:
        for j in range(0, narr.shape[0]):
          mline = narr[j,:]
          # add on posterior
          mline = np.append(mline, np.exp(mline[0]))
          
          strmline = " ".join(str(x) for x in mline) + '\n'
          cf.write(strmline)
      cf.close()
          
      # read in as common object
      peparser = bppu.PEOutputParser('common')
      cf = open(comfile, 'r')
      commonResultsObj = peparser.parse(cf)
      cf.close()
          
      # remove temporary file
      os.remove(comfile)
          
      # create posterior class
      pos = bppu.Posterior( commonResultsObj, SimInspiralTableEntry=None, \
                          votfile=None )
      
      # convert iota back to cos(iota)
      # create 1D posterior class of cos(iota) values
      cipos = bppu.PosteriorOneDPDF('cosiota', np.cos(pos['iota'].samples))
      
      # add it back to posterior
      pos.append(cipos)
      
      # remove iota samples
      pos.pop('iota')
      
      # get from info about the posteriors
      mP, mL = pos.maxL # maximum likelihood/posterior point
      sd = pos.stdevs #  standard deviations of posteriors 
      mv = pos.means # means of posteriors
      md = pos.medians # medians of posteriors
      cc, hn = corrcoef(pos) # correlation coefficients of posteriors
      
      poslist.append(pos) # append posterior object to list
      
      max_pos.append(mP)
      maxLvals.append(mL)
      stddevvals.append(sd)
      meanvals.append(mv)
      medvals.append(md)
      corrcoefs.append(cc)
      corrcoefsheader.append(hn)
    
    # if error occured in reading in MCMC chains then move to next pulsar
    if erridx:
      continue
    
    # set whether to attempt to output injection parameter on plot
    parinj = None
    if swinj or hwinj:
      parinj = par
      
    # output the MAIN posterior plots
    # h0
    bounds = [0, float("inf")]
    ul = 0.95
    histbins = 30
    h0Fig, ulvals = pppu.plot_posterior_hist( poslist, 'h0', ifosNew, \
                                              bounds, histbins, \
                                              ul, overplot=True, \
                                              parfile=parinj )
    h0figname = output_fig(h0Fig[0], puldir, 'h0post', ftypes)
    
    # convert h0 uls into ellipticity and spin-down ratio
    ell = []
    sdrat = []
    h0last = ulvals[-1]
    ullist.append(ulvals)
    
    limittable = []
    limittable.append( \
"""
<table>
<tr>
  <td colspan="4" style="text-align: center;">h<sub>0</sub> spin-down = %s</td>
</tr>
<tr>
  <th>&nbsp;</th>
  <th>h<sub>0</sub><sup>95%%</sup></th>
  <th>&#949;</th>
  <th>ratio</th>
</tr>
""" % exp_str(sdlim) )

    for i, ifo in enumerate(ifosNew):
      ell.append(pppu.h0_to_ellipticity(ulvals[i], f0, dist))
      sdrat.append(ulvals[i]/sdlim)
    
      # output values to main page table
      htmltext.append( \
"""
<td class="leftborder">%s</td>
<td>%s</td>
<td>%.2f</td>
""" % (exp_str(ulvals[i]), exp_str(ell[i]), sdrat[i]) )

      # output values to the main LaTeX table
      latextext.append( \
"""\
& $%s$ & $%s$ & $%.2f$ \
""" % (exp_latex_str(ulvals[i]), exp_latex_str(ell[i]), sdrat[i]))

      # output values to pulsar page table
      limittable.append( \
"""
<tr>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
  <td class="%s">%.2f</td>
</tr>
""" % (ifo, ifo, ifo, exp_str(ulvals[i]), ifo, exp_str(ell[i]), ifo, sdrat[i]) )

    htmltext.append('</tr>')
    latextext.append('\\\\')
    limittable.append('</table>')
    
    # phi0
    bounds = [0, 2*math.pi]
    phi0Fig, ulvals = pppu.plot_posterior_hist( poslist, 'phi0', ifosNew, \
                                        bounds, histbins, \
                                        0, overplot=True, parfile=parinj )
    phi0figname = output_fig(phi0Fig[0], puldir, 'phi0post', ftypes)
    
    # cos(iota)
    bounds = [-1, 1]
    ciFig, ulvals = pppu.plot_posterior_hist( poslist, 'cosiota', ifosNew, \
                                      bounds, histbins, \
                                      0, overplot=True, parfile=parinj )
    cifigname = output_fig(ciFig[0], puldir, 'cipost', ftypes)
    
    # psi
    bounds = [-math.pi/4, math.pi/4]
    psiFig, ulvals = pppu.plot_posterior_hist( poslist, 'psi', ifosNew, \
                                       bounds, histbins, \
                                       0, overplot=True, parfile=parinj )
    psifigname = output_fig(psiFig[0], puldir, 'psipost', ftypes)
   
    # get h0 vs cos(iota) 2D posterior histrogram (if single detector of for
    # joint posterior)
    #bounds = [[0, 3*h0last], [-1, 1]]
    h0ciFig = pppu.plot_posterior_hist2D([poslist[-1]], ['h0', 'cosiota'], \
[ifosNew[-1]], bounds=None, nbins=[30, 30], parfile=parinj)
    h0cifigname = output_fig(h0ciFig[0], puldir, 'h0cipost', ftypes)
    
    # get phi0 vs psi 2D posterior histogram
    phi0psiFig = pppu.plot_posterior_hist2D([poslist[-1]], ['phi0', 'psi'], \
[ifosNew[-1]], bounds=[[0, 2.*math.pi], [math.pi/4., -math.pi/4]], \
nbins=[30, 30], parfile=parinj)
    phi0psifigname = output_fig(phi0psiFig[0], puldir, 'phi0psipost', ftypes)
    
    # produce output table of posterior plots
    posttabletext = []
    posttabletext.append( \
"""
  <table>
    <tr>
      <td><a href="%s"><img class="posplot" src="%s"/></a></td>
      <td>&nbsp;</td>
      <td><a href="%s"><img class="posplot" src="%s"/></a></td>
      <td>&nbsp;</td>
    </tr>
    <tr>
      <td><a href="%s"><img class="posplot" src="%s"/></a></td>
      <td><a href="%s"><img class="rotated90clockwise" src="%s"/></a></td>
      <td><a href="%s"><img class="posplot" src="%s"/></a></td>
      <td><a href="%s"><img class="rotated90clockwise" src="%s"/></a></td>
    </tr>
  </table>
""" % (h0figname['png'], h0figname['png'], phi0figname['png'], \
phi0figname['png'], h0cifigname['png'], h0cifigname['png'], cifigname['png'], \
cifigname['png'], phi0psifigname['png'], phi0psifigname['png'], \
psifigname['png'], psifigname['png'] ))

    # get analysis statistics
    analysisstatstext = []
    analysisstatstext.append('<h2>Analysis statistics</h2>')

    analysisstatstext.append('<table>')
    for i, ifo in enumerate(ifos):
      # get the start time, end time and length from the heterodyned data file
      st = ((os.popen('head -1 %s' % Bkdata[i], "r")).readline()).split()[0]
      et = ((os.popen('tail -1 %s' % Bkdata[i], "r")).readline()).split()[0]
      lt = ((os.popen('wc -l %s' % Bkdata[i], "r")).readline()).split()[0]
      
      # duty cycle
      dc = 100.*(60.*float(lt))/(float(et)-float(st))
      
      analysisstatstext.append( \
"""
  <tr>
    <td style="text-align: center;"> 
      <table>
        <tr>
          <td>&nbsp;</td>
          <td class="%s">%s</td>
        </tr>
        <tr>
          <td>Start (GPS)</td>
          <td>%d</td>
        </tr>
        <tr>
          <td>End (GPS)</td>
          <td>%d</td>
        </tr>
        <tr>
          <td>Length (sec)</td>
          <td>%d</td>
        </tr>
        <tr>
          <td>Duty factor (%%)</td>
          <td>%.1f</td>
        </tr>
      </table>
    <td><a href="%s"><img class="dataplot" src="%s"/></a></td>
  </tr>
""" % (ifo, ifo,int(float(st)), int(float(et)), int(float(et)-float(st)), dc, \
Bkfigname[i]['png'], Bkfigname[i]['png']) )
  
      if plotpsds and plotfscan:
        analysisstatstext.append( \
"""
  <tr>
    <td><a href="%s"><img class="asdplot" src="%s"/></a></td>
    <td><a href="%s"><img class="dataplot" src="%s"/></a></td>
  </tr>
""" % (psdfigname[i]['png'], psdfigname[i]['png'], fscanfigname[i]['png'], \
fscanfigname[i]['png']) )
    analysisstatstext.append('</table>')
    
    # output MCMC chains and statistics
    mcmctabletext = []
    mcmctabletext.append('<h2>MCMC chains</h2>')
    mcmctabletext.append( \
"""<table>
  <tr>
    <th colspan=4>MCMC chain information</th>
  </tr>
  <tr>
    <th>&nbsp;</th>
    <th>no. of chains</th>
    <th>chain length</th>
    <th>effective sample size</th>
  </tr>
""" )
    for i, ifo in enumerate(ifosNew):
      mcmctabletext.append( \
"""
  <tr>
    <td class="%s">%s</td>
    <td>%d</td>
    <td>%d</td>
    <td>%d</td>
  </tr>
""" % (ifo, ifo, len(mcmcdirs), chainlens[i], nefflist[i]) )
    mcmctabletext.append('</table>')
    
    mcmctabletext.append('<table>')
    parlist = poslist[0].names
    for par in parlist:
      if par != 'post' and par != 'logl':
        chainfig = pppu.plot_posterior_chain(poslist, par, ifosNew, mcmcgr, \
                                             withhist=30)
        if chainfig:
          figname = output_fig(chainfig, puldir, par+'_mcmcchain', ftypes)
          mcmctabletext.append( \
"""
<tr>
<td><a href="%s"><img class="chainplot" src="%s"/></a></td>
</tr>
""" % (figname['png'], figname['png']) )
        else:
          break
    
    mcmctabletext.append('</table>')
    
    if chainfig==None:
      # if no MCMC chain figures were made (e.g. because gridspec wasn't
      # available) the remove the last table entries
      del mcmctabletext[-2:]
    
    # output info on the parameters (std dev, mean, max posterior etc)
    poststatstext = []
    poststatstext.append( \
"""\
<h2>
<a href="javascript:toggle('toggleText');">Posterior statistics</a>
</h2>\
""" )
    poststatstext.append('<div id="toggleText" style="display: none;">')
    poststatstext.append('<table>')
    
    poststatstext.append( \
"""\
  <tr>
    <th>&nbsp;</th>
    <th>&nbsp;</th>
    <th>max. posterior</th>
    <th>mean</th>
    <th>median</th>
    <th>&sigma;</th>
  </tr>\
""" )
    
    for j, param in enumerate(corrcoefsheader[0].split()):
      try:
        pdisp = paramtextdisp[param.upper()]
      except:
        pdisp = param
      
      poststatstext.append('<tr>\n<td rowspan="%d">%s</td>' % (len(ifosNew), \
pdisp) )
      
      try:
        dispfunc = paramdisp2.__dict__[param.upper()]
      except:
        dispfunc = paramdisp2.__dict__['DEFAULT'] 

      for i, ifo in enumerate(ifosNew):
        poststatstext.append('<td class="%s">%s</td>' % (ifo, ifo) )

        poststatstext.append('<td>%s</td>' % dispfunc(str(maxLvals[i][param])))

        poststatstext.append('<td>%s</td>' % dispfunc(str(meanvals[i][param])))
       
        poststatstext.append('<td>%s</td>' % \
dispfunc(str(medvals[i][param])))
    
        poststatstext.append('<td>%s</td>' % \
dispfunc(str(stddevvals[i][param])))
        
        poststatstext.append('</tr>')
        
    poststatstext.append('</table>\n</div>')
    
    # output footer giving how the file was made
    pfootertext = []
    
    pfootertext.append( \
"""<div id="footer">
%s - %s <br>
Command lines used:<br>
%s<br>
%s<br>
</div>
""" % (__author__, now.strftime('%a %d %b %Y'), p_args, __version__))
    
    pfootertext.append('</body>\n')
    pfootertext.append('</html>')
    
    # put all the html bits together in the right order
    htmlptext.append(cattext(pheadertext))
    
    # put parameter table and upper limit table together in a table element
    htmlptext.append('<div class="wrapper">')
    htmlptext.append('<table>\n<tr>\n<td>')
    htmlptext.append(cattext(partext))
    htmlptext.append('</td>\n<td>')
    htmlptext.append(cattext(limittable))
    htmlptext.append('</tr>\n</table>\n</div>')
    
    # output the posterior plots in an element
    htmlptext.append('<div class="wrapper">')
    htmlptext.append(cattext(posttabletext))
    htmlptext.append('</div>')
    
    # output the analysis stats
    htmlptext.append('<div class="wrapper">')
    htmlptext.append(cattext(analysisstatstext))
    htmlptext.append('</div>')

    # output the MCMC chains
    htmlptext.append('<div class="wrapper">')
    htmlptext.append(cattext(mcmctabletext))
    htmlptext.append('</div>')
    
    # output the posterior stats
    htmlptext.append('<div class="wrapper">')
    htmlptext.append(cattext(poststatstext))
    htmlptext.append('</div>')
    
    # output footer
    htmlptext.append(cattext(pfootertext))
    
    # add put it all together
    htmlpout.write(cattext(htmlptext))
    htmlpout.close()

  htmltext.append('</table>')

  htmltext.append( \
"""
<p>
  <sup>&dagger;</sup> The spin-down limit was calculated using a spin-down
  <a \
href="http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=\
normal&highlight=age_i#age_i">corrected</a> for proper motion effects.
</p>""" )
  
  # output footer giving how the file was made
  htmltext.append( \
"""<div id="footer">
%s - %s <br>
Command lines used:<br>
%s<br>
%s<br>
</div>
""" % (__author__, now.strftime('%a %d %b %Y'), p_args, __version__))
 
  htmltext.append('</body>\n') 
  htmltext.append('</html>\n')
  htmlout.write(cattext(htmltext))
  htmlout.close()
  
  latextext.append('\\enddata')
  if iscorlatex:
    latextext.append("\\tablenotetext{\dagger}{The pulsar's spin-down is \
corrected for proper motion effects}")
  latextext.append('\\end{deluxetable}')
  latextext.append('\clearpage')
  latextext.append('\\end{landscape}')
  latextext.append('\\end{document}')
  latexout.write(cattext(latextext))
  latexout.close()
  
  # compile the LaTeX table
  sp.call(['latex', latexpage])
  sp.call(['latex', latexpage])
  sp.call(['dvips', re.sub('tex', 'dvi', latexpage)])
  sp.call(['ps2pdf', re.sub('tex', 'ps', latexpage)])
  
