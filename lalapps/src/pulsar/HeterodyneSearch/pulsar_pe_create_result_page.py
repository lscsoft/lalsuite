#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       pulsar_mcmc_create_results_page.py
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

#standard library imports
import sys
import os
import re
import shutil
import math
import datetime
import shelve
import struct
import gzip
from optparse import OptionParser
import subprocess as sp

#related third party imports
import numpy as np

# need to set Agg here for use of stuff in pulsarpputils
import matplotlib
matplotlib.use("Agg")

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
  text-align: left;
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
  border-collapse: collapse;
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


# return the ATNF information on a given pulsar
def get_atnf_info(par):
  try:
    # set possible pulsar names to try consecutively, also then try each with the wildcard * (this is currently necessary for J2022+17)
    trynames = [par['PSRJ'], par['PSRB'], par['PSR'], par['NAME'], par['PSRJ'], par['PSRB'], par['PSR'], par['NAME']]

    for i, psrname in enumerate(trynames):
      badurl = False

      if psrname != None:
        psrname = re.sub('\+', '%2B', psrname) # switch + for unicode character

        if i+1 > len(trynames)/2:
          psrname = psrname + '*'

        atnfurl = \
'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=1.53&Dist=Dist&Assoc=\
Assoc&Age_i=Age_i&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&\
sort_attr=&sort_order=asc&condition=&pulsar_names=' + psrname + \
'&ephemeris=selected&submit_ephemeris=\
Get+Ephemeris&coords_unit=raj%2Fdecj&radius=&coords_1=\
&coords_2=&style=Long+with+last+digit+error&no_value=*&fsize=3&x_axis=&x_scale=\
linear&y_axis=&y_scale=linear&state=query'

        soup = None
        pdat = None

        soup = bs(urllib2.urlopen(atnfurl).read())
        pdat = soup.pre # data exists in the pre html environment

        for line in pdat:
          vals = line.split('\n') # split at any new lines

          for row in vals:
            if 'WARNING' in row or 'not in catalogue' in row:
              badurl = True
              break

          if badurl:
            break

        if badurl:
          continue
        else:
          return (pdat, psrname)

    if badurl:
      return None
  except:
    return None

# list of parameters to display (in this order)
paramdisplist = ['RAJ', 'DECJ', 'F0', 'F1', 'F2', 'PEPOCH', 'X', 'E', \
'EPS1', 'EPS2', 'OM', 'T0', 'TASC', 'PB']

# html text to display for different parameter names
paramtextdisp = {'RAJ': '&alpha;', 'DECJ': '&delta;', \
                 'RA': '&alpha;', 'DEC': '&delta;', \
                 'F0': 'f<sub>0</sub> (Hz)', 'F1': 'f<sub>1</sub> (Hz/s)', \
                 'F2': 'f<sub>2</sub> (Hz/s<sup>2</sup>)', \
                 'F3': 'f<sub>3</sub> (Hz/s<sup>3</sup>)', \
                 'F4': 'f<sub>4</sub> (Hz/s<sup>4</sup>)', \
                 'F5': 'f<sub>5</sub> (Hz/s<sup>5</sup>)', \
                 'PEPOCH': 'epoch (MJD)', 'A1': 'a sin<it>i</i> (lt s)', \
                 'E': '<it>e</it>', 'EPS1': '&epsilon;<sub>1</sub>', \
                 'EPS2': '&epsilon;<sub>2</sub>', \
                 'EPS1': '&epsilon;<sub>1</sub>', \
                 'T0': 'T<sub>0</sub> (MJD)', \
                 'TASC': 'T<sub>asc</sub> (MJD)', \
                 'OM': '&omega;<sub>0</sub> (deg)',
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
  def A1(f): return dec_or_exp(f)
  def E(f): return dec_or_exp(f)
  def EPS1(f): return exp_str(float(f), 2)
  def EPS2(f): return exp_str(float(f), 2)
  def T0(f): return dec_or_exp(f)
  def TASC(f): return dec_or_exp(f)
  def OM(f): return dec_or_exp(f)
  def PB(f): return dec_or_exp(f)
  def H0(f): return exp_str(float(f), 2)
  def COSIOTA(f): return '%.2f' % float(f)
  def PHI0(f): return '%.2f' % float(f)
  def PSI(f): return '%.2f' % float(f)
  def PMRA(f): return exp_str(float(f), 2)
  def PMDC(f): return exp_str(float(f), 2)
  def DEFAULT(f): return '%.2f' % float(f)


# function will return a string with the number (input as a string) to two decimal
# places if it is greater than 0.01, or the number in exponent for to one decimal
# place it it is smaller
def dec_or_exp(f):
  if float(f) > 0.01:
    return '%.2f' % float(f)
  else:
    return exp_str(float(f), 1)


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
  description = \
"""This script is for creating a results output page and plots for the known
   pulsar analysis. It can use either the inputs from the old MCMC code
   lalapps_pulsar_parameter_estimation, or the nested sampling code
   lalapps_pulsar_parameter_estimation_nested, as requested."""

  epilog = \
"""
An example of usage for a case when two MCMC runs have been
performed for two interferometers (H1 and L1):
  %s --ifo H1 --Bkfiles /home/me/hetdata/H1 --ifo L1 --Bkfiles
/home/me/hetdata/L1 --parfile /home/me/pardir/pulsar.par --priorfile priors.txt
--histbins 50 --mcmc --mcmcdirs /home/me/MCMCchains/pdfs1
--mcmcdirs /home/me/MCMCchains/pdfs2 --outpath /home/me/public_html/results
""" % os.path.basename(sys.argv[0])
  usage = "Usage: %prog [options]"

  parser = OptionParser( usage = usage, description = description,
                         version = __version__, epilog = epilog )

  parser.add_option("-o", "--outpath", dest="outpath", help="The path for "
                    "the analysis output (a sub-directory based on the pulsar "
                    "name will be created here that contains the pulsar "
                    "information)", metavar="DIR")

  """
   data will be read in in the following format:
     multiple MCMC directories from each IFO would be entered as follows:
     --ifo H1 --Bkfiles /home/user/finehetsH1
     --ifo L1 --Bkfiles /home/user/finehetsL1
     --mcmcdirs /home/user/pdfs1
     --mcmcdirs /home/user/pdfs2
     --mcmcdirs /home/user/pdfs3
  """
  parser.add_option("-M", "--mcmc", dest="usemcmc",
                    help="Input files will be from the MCMC code.",
                    action="store_true", default=False)

  parser.add_option("-m", "--mcmcdirs", dest="mcmcdirs",\
                    action="append", help="If using the MCMC inputs supply an \
MCMC directory containing chains for each IFO.")

  parser.add_option("-N", "--nested", dest="usenested",
                    help="Input files will be samples from the nested \
sampling code and pre-converted into posterior samples through lalapps_nest2pos.",
                    action="store_true", default=False)

  parser.add_option("-f", "--nestedfiles", dest="nestedfiles",
                    action="append", help="If using nested sampling inputs \
include one posterior sample file for each IFO.")

  # get pulsar .par file
  parser.add_option("-p", "--parfile", dest="parfile", help="An "
                    "individual, TEMPO-style pulsar parameter file, used in "
                    "the analysis. [required]", metavar="PNAME.par",
                    default=None)

  # get heterodyned data files used for analysis
  parser.add_option("-b", "--Bkfiles", dest="Bkfiles", action="append",
                    help="A directory containing the heterodyned data files "
                    "for a given detector [required]", default=None)

  # get the prior file used for the analysis
  parser.add_option("-r", "--priordir", dest="priorfile", help="A directory "
                    "containing the 2D prior files, on h0 and cos(iota), used "
                    "in the analysis.", default=None)

  # get list of the detectors used (same order as Bk files)
  parser.add_option("-i", "--ifo", dest="ifos", help="The individual "
                    "interferometers from which analysis output, and data "
                    "files have been supplied. [required]", action="append",
                    default=None)

  # get number of bins for histogramming (default = 50)
  parser.add_option("-n", "--histbins", dest="histbins", help="The number of " \
                    "bins for histrogram plots. [default = 50]", default=50)

  # say if you want to output eps figure (default to False)
  parser.add_option("-e", "--epsout", dest="epsout", help="Set if wanting to " \
                    "output eps versions of figures", action="store_true",
                    default=False)

  # say if the pulsar is a software injections (default to False)
  parser.add_option("-s", "--sw-inj", dest="swinj", help="Set if the " \
                    "pulsars are software injections", action="store_true", \
                    default=False)

  # say if the pulsar is a hardware injections (default to False)
  parser.add_option("-w", "--hw-inj", dest="hwinj", help="Set if the " \
                    "pulsars are hardware injections", action="store_true", \
                    default=False)

  # get confidence interval
  parser.add_option("-c", "--confidence-interval", dest="ci", help="Set the \
 % confidence interval to calculate [default: %default %].", type="float",
                    default=95.)

  p_args=''
  for arg in sys.argv:
    p_args += arg + ' '

  print >> sys.stderr, p_args

  # parse input options
  (opts, args) = parser.parse_args()

  # check that output path has been given
  if not opts.__dict__['outpath']:
    print >> sys.stderr, "Must specify an output path"
    parser.print_help()
    sys.exit(0)
  else:
    outpath = opts.outpath

  usemcmc = opts.usemcmc
  usenested = opts.usenested

  # check that some data has been given
  if usemcmc:
    if not opts.__dict__['mcmcdirs']:
      print >> sys.stderr, "Must specify MCMC chain directories"
      parser.print_help()
      sys.exit(0)
    else:
      # split MCMC directories
      mcmcdirs = opts.mcmcdirs
  elif usenested:
    if not opts.__dict__['nestedfiles']:
      print >> sys.stderr, "Must specify nested sampling files."
      parser.print_help()
      sys.exit(0)
    else:
      nestedfiles = opts.nestedfiles
  else:
    print >> sys.stderr, "Must specify using either the MCMC or nested sampling inputs."
    parser.print_help()
    sys.exit(0)

  # check that parfile has been given
  if not opts.__dict__['parfile']:
    print >> sys.stderr, "Must specify a pulsar TEMPO .par file"
    parser.print_help()
    sys.exit(0)

  if not opts.__dict__['ifos']:
    print >> sys.stderr, "Must specify the interferometers analysed"
    parser.print_help()
    sys.exit(0)
  else:
    ifos = opts.ifos

  if opts.__dict__['Bkfiles'] is None:
    Bkfiles = []
  else:
    Bkfiles = opts.Bkfiles

  if opts.epsout:
    ftypes = ['png', 'eps']
  else:
    ftypes = ['png']

  priordir = None
  if opts.priorfile:
    if os.path.isdir(opts.priorfile):
      priordir = opts.priorfile

  swinj = opts.swinj
  hwinj = opts.hwinj

  # create output directory if it doesn't exist
  if not os.path.isdir(outpath):
    try:
      os.mkdir(outpath)
    except:
      print >> sys.stderr, "Cannot create output directory %s" % outpath
      sys.exit(0)

  # check that number of ifos is the same as the number of data lists (for MCMC)
  nifos = len(ifos)
  ndata = len(Bkfiles)

  ifosNew = list(ifos)
  if usemcmc:
    if nifos > 1:
      ifosNew.append('Joint')
  if usenested:
    ifostmp = []

    # create shortened list of ifos removing the Joint one if specified
    if nifos > 1:
      for i in ifos:
        if 'Joint' in i:
          continue
        else:
          ifostmp.append(i)
      ifos = ifostmp
      nifos = len(ifos)

  # if we only want a Joint output then set this flag to true - this will prevent the code requiring Bk files
  jointonly = False
  if nifos == 1 and ifos[0] == 'Joint':
    jointonly = True

  # check that number of ifos is the same as the number of data lists (for MCMC)
  if not jointonly: # if we only want a Joint output don't perform this check
    if nifos != ndata:
      print >> sys.stderr, "Number of IFOs and data lists are not equal"
      sys.exit(0)

  # check if parfile is a single file or a directory
  if os.path.isfile(opts.parfile):
    parfile = opts.parfile
  else:
    print >> sys.stderr, "No par file specified!"
    sys.exit(0)

  # get time/date for file creation
  now = datetime.datetime.now()

  # list of amplitude spectral densities, upper limits and spin-down ratios
  # for the IFOs
  asds = []
  sdlist = []

  # read data from par file
  try:
    par = pppu.psr_par(parfile)
  except:
    print >> sys.stderr, "Par file %s could not be opened!" % parfile
    sys.exit(0)

  pname = par['PSRJ']
  if not pname:
    print >> sys.stderr, "No PSRJ value in par file %s" % parfile

    pname = par['PSR']
    if not pname:
      print >> sys.stderr, "No PSR value in par file %s" % parfile
      sys.exit(0)

  # check if heterodyned data exists for this pulsar and each detector
  if not jointonly:
    Bkdata = []
    for i, ifo in enumerate(ifos):
      Bkdata.append(os.path.join(Bkfiles[i], 'finehet_' + pname + '_' + ifo))

      # check files exist if not then skip the pulsar
      if not os.path.isfile(Bkdata[i]) and not os.path.isfile(Bkdata[i]+'.gz'):
        print >> sys.stderr, "No heterodyne file %s" % Bkdata[i]
        sys.exit(0)

  # check that MCMC chains exist for this pulsar and each detector (including
  # joint)
  if usemcmc:
    for ifo in ifosNew:
      for chaindir in mcmcdirs:
        cfile = os.path.join(chaindir, 'MCMCchain_' + pname + '_' + ifo)

        if not os.path.isfile(cfile):
          print >> sys.stderr, "No MCMC file %s" % cfile
          sys.exit(0)
  if usenested:
    if len(ifosNew) != len(nestedfiles):
      print >> sys.stderr, "Number of nested sampling file lists must be equal to number of IFOs."
      sys.exit(-1)

    # check files exist
    for nfile in nestedfiles:
      if not os.path.isfile(nfile):
        print >> sys.stderr, "No nested sampling file %s" % nfile
        sys.exit(0)

  # check required parameters
  f0 = par['F0']
  if not f0:
    print >> sys.stderr, "No F0 value in par file %s" % parfile
    sys.exit(0)

  f1 = par['F1']
  if not f1:
    print >> sys.stderr, "No F1 value in par file %s, setting to zero" % parfile
    f1 = 0.

  print 'Results for pulsar ' + pname

  # create output directory for pulsar
  puldir = os.path.join(outpath, pname)
  if not os.path.isdir(puldir):
    os.mkdir(puldir)

  # create shelf database to store info for later parsing
  try:
    psrshelf = shelve.open(os.path.join(puldir, pname+'.db'))
  except:
    print >> sys.stderr, "Can't open shelf database!"
    sys.exit(1)

  psrshelf['par'] = par
  psrshelf['ifos'] = ifosNew

  # create CSS
  cssname = 'psr.css'
  try:
    cssfile = os.path.join(puldir, cssname)
    css = open(cssfile, 'w')
  except:
    print >> sys.stderr, "Cannot open CSS file %s" % cssfile
    sys.exit(1)

  css.write(csstext);
  css.close()

  # attempt to get pulsar distance, proper motion corrected age and any
  # association (e.g. GC) from the ATNF catalogue
  ages = None
  dists = None
  assoc = None
  notinatnf = False # flag to check if pulsar is in the ATNF catalogue or not
  age_i = None
  intrinsicsd = False
  agebasedsd = False

  if not swinj and not hwinj:
    atnfinfo = get_atnf_info(par)

    if atnfinfo != None:
      pdat = atnfinfo[0]
      pnameurl = atnfinfo[1]

      for line in pdat:
        vals = line.split('\n') # split at any new lines

        for row in vals:
          if 'DIST' in row:
            dists = row.split() # split row at whitespace

          if 'AGE_I' in row:
            ages = row.split()

          if 'ASSOC' in row:
            assoc = row.split()
    else:
      notinatnf = True
      atnfurl = None
      print >> sys.stderr, 'Problem accessing ATNF for %s!' % pname

  if ages:
    if len(ages) > 1:
      try:
        age_i = float(ages[1])
        intrinsicsd = True # using intrinsic spin-down
      except:
        print >> sys.stderr, "%s: Age is not a number!" % pname
        age_i = None
        f1sd = None

      f1sd = -f0/(2*age_i*365.25*86400)
    else:
      age_i = None
      f1sd = f1
  # check if pulsar is a GC pulsars for which no proper motion corrected age
  # is known and for these cases set a conservative upper limit assuming
  # tau = f0/(2*f1) = 10^9 years
  elif assoc:
    if len(assoc) > 1:
      if 'GC' in assoc[1]:
        f1sd = -f0/(2*1e9*365.25*86400) # set f1 for spin-down limit calc
        agebasedsd = True # using a conservative age-based spin-down limit
      else:
        f1sd = f1 # some Fermi pulsars have an assoc, but are no in GCs
                  # Note: also currently (21/01/13) some Fermi pulsars
                  # wrongly are given a spin-up in the ATNF catalogue!
  else:
    f1sd = f1

  psrshelf['f0rot'] = f0
  psrshelf['f0gw'] = f0*2.
  psrshelf['f1sd'] = f1sd
  psrshelf['f1rot'] = f1sd
  psrshelf['f1gw'] = f1sd*2.
  psrshelf['age_i'] = age_i
  psrshelf['intrinsicsd'] = intrinsicsd
  psrshelf['agebasedsd'] = agebasedsd

  # first try getting a distance from a par file
  dist = par['DIST']

  # if pulsar is not in atnf and there's not distance in the par file then skip it
  if notinatnf and not dist:
    print >> sys.stderr, "No distance available for %s" % pname

  # if not available look in download value
  if not notinatnf and not dist:
    # see if it was in the ATNF catalogue
    if dists:
      if len(dists) > 1:
        try:
          dist = float(dists[1])
        except:
          print >> sys.stderr, "%s: Distance not a number!" % pname
          dist = None
      else:
        print >> sys.stderr, "%s: Distance not a number!" % pname
        dist = None

  psrshelf['dist'] = dist

  # create html page
  htmlppage = os.path.join(puldir, 'index.html')
  try:
    htmlpout = open(htmlppage, 'w')
  except:
    print >> sys.stderr, 'Could not open html page %s' % htmlppage

  htmlptext = None
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

  <link rel="stylesheet" type="text/css" href="%s"/>

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
  pheadertext = None
  pheadertext = []

  # set title prefix depending is injection or not
  psrnameprefix = 'PSR'
  if swinj:
    psrnameprefix = 'SWINJ'
  elif hwinj:
    psrnameprefix = 'HWINJ'

  if not swinj and not hwinj and not notinatnf:
    atnfurl = 'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=1.47&startUserDefined=true&pulsar_names=' + pnameurl + '&ephemeris=long&submit_ephemeris=Get+Ephemeris&state=query'
    pheadertext.append('<h1><a href="%s">%s %s</a></h1>\n' %
(atnfurl, psrnameprefix, pname))
  else:
    pheadertext.append('<h1>%s %s</h1>\n' % (psrnameprefix, pname))

  # copy par file to pulsar directory
  shutil.copy(parfile, os.path.join(puldir, pname + '.par'))

  # print out par file info
  partext = None
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

  psrshelf['ra'] = par['RAJ']
  psrshelf['dec'] = par['DECJ']

  # get spin-down upper limit
  if not swinj and not hwinj and f1sd and dist:
    sdlim = pppu.spin_down_limit(f0, f1sd, dist)
  else:
    sdlim = None

  psrshelf['sdlim'] = sdlim

  # get time series and PSD plots
  if not jointonly:
    Bkdata = []
    plotpsds = True
    plotfscan = True

    for i, ifo in enumerate(ifos):
      Bkdata.append(Bkfiles[i] + '/finehet_' + pname + '_' + ifo)
      # check file exists
      if not os.path.isfile(Bkdata[i]):
        Bkdata[i] = Bkdata[i]+'.gz' # try gzipped file
        if not os.path.isfile(Bkdata[i]):
          print >> sys.stderr, "Error... could not find Bk data file %s" % Bkdata[i]

    asdtime = 14400 # set time over which to produce the asds

    Bkfigs, psdfigs, fscanfigs, asdlist = pppu.plot_Bks_ASDs( Bkdata, ifos, \
asdtime, plotpsds=plotpsds, plotfscan=plotfscan, removeoutlier=50 )

    if asdlist:
      psrshelf['ASD'] = dict(zip(ifos, asdlist)) # convert into dictionary

    # output plots of time series and psd
    Bkfigname = None
    psdfigname = None
    fscanfigname = None

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
  evidence = []
  confidenceinterval = []

  if swinj or hwinj:
    confidenceregion = []

  for i, ifo in enumerate(ifosNew):
    if usemcmc:
      # read in MCMC chains
      chainfiles = []

      for chaindir in mcmcdirs:
        chainfiles.append(os.path.join(chaindir, 'MCMCchain_'+pname+'_'+ifo))

      pos, neffs, grr, cl = pppu.pulsar_mcmc_to_posterior(chainfiles)

      if pos == None or neffs == None or grr == None or cl == None:
        outerrfile = os.path.join(outpath, 'mcmc_chain_errors.txt')
        errfile = open(outerrfile, 'a') # append to file
        errfile.write('MCMC files for %s and %s could not be read in\n' % (pname, ifo))
        errfile.close()
        print >> sys.stderr, "Error in MCMC chain reading!"

        try:
          shutil.rmtree(puldir)
        except:
          os.system('rm -rf %s' % puldir)

        sys.exit(0)

      nefflist.append(math.fsum(neffs))
      chainlens.append(np.mean(cl))
      mcmcgr.append(grr)

    if usenested:
      pos, ev = pppu.pulsar_nest_to_posterior(nestedfiles[i])
      evidence.append(ev)
      mcmcgr = None
    else:
      evidence.append(None)

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

    # get confidence interval
    cidict = {} # dictionary of confidence intervals
    if swinj or hwinj:
      ciregion = {} # dictionary of confidence interval containing injection value
    pnames = pos.names
    # remove 'logL', 'deltalogl', 'deltalogL', 'logw' or 'logPrior'
    for pn in pnames:
      pl = pn.lower()
      if 'logl' in pl or 'deltalogl' in pl or 'logw' in pl or 'logprior' in pl:
        continue
      else:
        # DONT NEED TO USE KDTREE STUFF FOR THIS (BUT KEEP HERE FOR REFERENCE)
        #leaves, intArea, injInfo = bppu.kdtree_bin2step(pos, [pn], [opts.ci/100.], samples_per_bin=64)

        #lowbound = np.inf
        #highbound = -np.inf

        #for j in range(len(leaves)):
        #  if leaves[j][4] > opts.ci/100.:
        #    break # found edge of bound
        #  else:
        #    lowbound = np.amin([leaves[j][0][0][0], leaves[j][0][1][0], lowbound])
        #    highbound = np.amax([leaves[j][0][0][0], leaves[j][0][1][0], highbound])

        psamples = np.squeeze(pos[pn].samples).tolist()
        psamples.sort() # sort samples into ascending order
        lowbound = min(psamples)
        highbound = max(psamples)
        cregion = highbound - lowbound
        lsam = len(psamples)
        cllen = int(lsam*opts.ci/100.)
        for j in range(lsam-cllen):
          if psamples[j+cllen]-psamples[j] < cregion:
            lowbound = psamples[j]
            highbound = psamples[j+cllen]
            cregion = highbound - lowbound

        # get confidence region containing the injection
        if swinj or hwinj:
          # loop over difference confidence intervals until finding the smallest one
          # that contains the injection
          cifound = False
          for ci in range(1,101):
            lowr = min(psamples)
            highr = max(psamples)
            cregion = highr - lowr
            cllen = int(lsam*ci/100.)
            for j in range(lsam-cllen):
              if psamples[j+cllen]-psamples[j] < cregion:
                lowr = psamples[j]
                highr = psamples[j+cllen]
                cregion = highr - lowr

            if par[pn.upper()] >= lowr and par[pn.upper()] <= highr:
              # found injection confidence region
              cifound = True
              ciregion[pn] = ci
              break

          if not cifound:
            ciregion[pn] = 101. # injection is outside of posterior!

        cidict[pn] = [lowbound, highbound]

    confidenceinterval.append(cidict)

    if swinj or hwinj:
      confidenceregion.append(ciregion)

  # set whether to attempt to output injection parameter on plot
  parinj = None
  if swinj or hwinj:
    parinj = par

  # output the MAIN posterior plots
  # h0
  bounds = [0, float("inf")]
  ul = 0.95
  if not opts.__dict__['histbins']:
    histbins = 30 # default number of histogram bins
  else:
    histbins = opts.histbins

  h0Fig, ulvals = pppu.plot_posterior_hist( poslist, 'h0', ifosNew, \
                                            bounds, histbins, \
                                            ul, overplot=True, \
                                            parfile=parinj )
  h0figname = output_fig(h0Fig[0], puldir, 'h0post', ftypes)

  # convert h0 uls into ellipticity and spin-down ratio
  ell = []
  sdrat = []
  sdpowrat = []
  q22 = []
  h0last = ulvals[-1]

  psrshelf['h95'] = dict(zip(ifosNew, ulvals)) # convert to dictionary

  limittable = []

  sdtext = ''
  sdouttext = ''
  if intrinsicsd:
    sdtext = '<sup>&dagger;</sup>' # dagger
    sdouttext = '<p>%sThe spin-down limit was calculated using a spin-down corrected for proper motion effects.</p>' % sdtext
  elif agebasedsd:
    sdtext = '<sup>&Dagger;</sup>' # double dagger
    sdouttext = '<p>%sThe spin-down limit was calculated using a characteristic age of 10<sup>9</sup> years</p>' % sdtext

  if sdlim == None:
    sdlimtxt = ''
  else:

    sdlimtxt = \
"""\
<tr>
  <td colspan="6" style="text-align: center;">h<sub>0</sub> spin-down = %s%s</td>
</tr>\
""" % (exp_str(sdlim), sdtext)

  limittable.append( \
"""
<table>
%s
<tr>
  <th>&nbsp;</th>
  <th>h<sub>0</sub><sup>95%%</sup></th>
  <th>&#949;</th>
  <th>Q<sub>22</sub></th>
  <th>ratio</th>
  <th>log(evidence ratio)</th>
</tr>
""" % sdlimtxt )

  for i, ifo in enumerate(ifosNew):
    if not swinj and not hwinj and dist:
      ell.append(pppu.h0_to_ellipticity(ulvals[i], f0, dist))
      ellhtmlstr = '%s' % exp_str(ell[i])
      sdrat.append(ulvals[i]/sdlim);
      sdratstr = '%.2f' % sdrat[i]
      q22.append(pppu.h0_to_quadrupole(ulvals[i], f0, dist))
      q22htmlstr = '%s' % exp_str(q22[i])
      sdpowrat.append(100.*sdrat[i]**2)
    else:
      ell.append(None)
      q22.append(None)
      sdrat.append(None)
      sdpowrat.append(None)
      ellhtmlstr = '*'
      q22htmlstr = '*'
      sdratstr = '*'

    if usenested:
      if evidence[i] > 100. or evidence[i] < -100.:
        evidencestr = '%s' % exp_str(evidence[i])
      else:
        evidencestr = '%.2f' % evidence[i]
    else:
      evidencestr = '*'

    # output values to pulsar page table
    limittable.append( \
"""
<tr>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
  <td class="%s">%s</td>
</tr>
""" % (ifo, ifo, ifo, exp_str(ulvals[i]), ifo, ellhtmlstr, ifo, q22htmlstr, ifo, sdratstr, ifo, evidencestr) )

  limittable.append('</table>')

  psrshelf['ell'] = dict(zip(ifosNew, ell))
  psrshelf['q22'] = dict(zip(ifosNew, q22))
  psrshelf['sdrat'] = dict(zip(ifosNew, sdrat))
  psrshelf['sdpowrat'] = dict(zip(ifosNew, sdpowrat))
  psrshelf['evidence'] = dict(zip(ifosNew, evidence))
  psrshelf['confidenceinterval'] = dict(zip(ifosNew, confidenceinterval))
  psrshelf['ci'] = opts.ci # the percentage confidence interval

  if swinj or hwinj:
    psrshelf['confidenceregion'] = dict(zip(ifosNew, confidenceregion))

  # phi0
  bounds = [0, math.pi]
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
  bounds = [0., math.pi/2.]
  psiFig, ulvals = pppu.plot_posterior_hist( poslist, 'psi', ifosNew, \
                                       bounds, histbins, \
                                       0, overplot=True, parfile=parinj )
  psifigname = output_fig(psiFig[0], puldir, 'psipost', ftypes)

  # get h0 vs cos(iota) 2D posterior histrogram (if single detector of for
  # joint posterior)
  # get bounds from h0 and cos(iota) plots
  ax = h0Fig[0].get_axes()
  gh0 = ax[-1].axis()
  ax = ciFig[0].get_axes()
  gci = ax[-1].axis()
  h0cibounds = [[gh0[0], gh0[1]], [gci[1], gci[0]]]

  h0ciFig = pppu.plot_posterior_hist2D(poslist, ['h0', 'cosiota'], \
ifosNew, bounds=h0cibounds, nbins=[30, 30], parfile=parinj, overplot=True)
  h0cifigname = output_fig(h0ciFig[0], puldir, 'h0cipost', ftypes)

  # get phi0 vs psi 2D posterior histogram
  ax = phi0Fig[0].get_axes()
  gphi0 = ax[-1].axis()
  ax = psiFig[0].get_axes()
  gpsi = ax[-1].axis()
  phi0psibounds = [[gphi0[0], gphi0[1]], [gpsi[1], gpsi[0]]]

  phi0psiFig = pppu.plot_posterior_hist2D(poslist, ['phi0', 'psi'], \
ifosNew, bounds=phi0psibounds, nbins=[30, 30], parfile=parinj, overplot=True)
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

  # get prior file if available
  h0prior = None
  priortabletext = []
  if priordir:
    priorfile = os.path.join(priordir, 'h0cipost_'+pname+'.bin')

    if os.path.isfile(priorfile):
      figprior = pppu.plot_2Dhist_from_file(priorfile, 'h0', 'cosiota')

      # get the previous h0 upper limit from that prior file
      h0prior = pppu.h0ul_from_prior_file(priorfile, ulval=0.95)
      print h0prior
    else:
      figprior = h0prior = None

    if not figprior or not h0prior:
      print >> sys.stderr, "Could not create prior file figures"
      priorfile = None
      h0prior = None
    else:
      # output plots
      priorh0cifigname = output_fig(figprior[0], puldir, 'h0ciprior', ftypes)
      priorh0figname = output_fig(figprior[1], puldir, 'h0prior', ftypes)
      priorcifigname = output_fig(figprior[2], puldir, 'ciprior', ftypes)

      priortabletext.append('<h2>Prior data</h2>')
      priortabletext.append('<p>h<sub>0</sub><sup>95%%</sup> prior upper limit: %s</p>' % exp_str(h0prior))
      priortabletext.append( \
"""
  <table>
    <tr>
      <td><a href="%s"><img class="posplot" src="%s"/></a></td>
      <td>&nbsp;</td>
    </tr>
    <tr>
      <td><a href="%s"><img class="posplot" src="%s"/></a></td>
      <td><a href="%s"><img class="rotated90clockwise" src="%s"/></a></td>
    </tr>
  </table>
""" % (priorh0figname['png'], priorh0figname['png'], priorh0cifigname['png'], \
priorh0cifigname['png'], priorcifigname['png'], priorcifigname['png']))

  psrshelf['h0prior'] = h0prior

  # get analysis statistics
  if not jointonly:
    analysisstatstext = []
    analysisstatstext.append('<h2>Analysis statistics</h2>')

    analysisstatstext.append('<table>')

    ulestimatetop = [] # list to contain upper limit estimates based on ASDs and data spans
    ulestimatebot = []
    lt = [] # length of the data for an IFO

    for i, ifo in enumerate(ifos):
      # get the start time, end time and length from the heterodyned data file
      st = et = None
      if '.gz' not in Bkdata[i]: # just use Popen to get data from heterodyne file
        st = (sp.Popen(['head', '-1', Bkdata[i]], stdout=sp.PIPE).communicate()[0]).split()[0]
        et = (sp.Popen(['tail', '-1', Bkdata[i]], stdout=sp.PIPE).communicate()[0]).split()[0]
        lt.append(float((sp.Popen(['wc', '-l', Bkdata[i]], stdout=sp.PIPE).communicate()[0]).split()[0])*60)
      else: # otherwise if gzipped we'll have to open the file
        bkd = np.loadtxt(Bkdata[i])
        st = bkd[0,0]
        et = bkd[-1,0]
        lt.append(float(len(bkd))*60.)

      # duty cycle
      dc = 100.*lt[i]/(float(et)-float(st))

      # get UL estimate based on h95 ~ 7-20 *sqrt(2) * ASD / sqrt(T) - the sqrt(2) is because
      # the ASD is calulated from a two-sided PSD
      if asdlist:
        ulestimatebot.append(7.*math.sqrt(2.)*np.median(asdlist[i])/math.sqrt(lt[i]))
        ulestimatetop.append((20./7.)*ulestimatebot[i])

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

    # get joint upper limit estimate
    if asdlist:
      if len(ifos) > 1:
        uli = []
        for i, asdv in enumerate(asdlist):
          uli.append(lt[i]/np.median(asdv)**2)

        ulesttmp = math.sqrt(2.)*math.sqrt(1./sum(uli))
        ulestimatebot.append(7.*ulesttmp)
        ulestimatetop.append(20.*ulesttmp)

    psrshelf['ulestimatebot'] = dict(zip(ifosNew, ulestimatebot))
    psrshelf['ulestimatetop'] = dict(zip(ifosNew, ulestimatetop))

  # output MCMC chains and statistics
  mcmctabletext = None
  mcmctabletext = []

  if usemcmc:
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
    <th>mean effective sample size</th>
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

    fnamepostfix = '_mcmcchain'

  if usenested:
    mcmctabletext.append('<h2>Posterior samples</h2>')
    fnamepostfix = '_postchain'

  mcmctabletext.append('<table>')
  parlist = None
  parlist = poslist[0].names
  for param in parlist:
    pl = param.lower()

    if 'post' in pl or 'logl' in pl or 'logw' in pl or 'logprior' in pl:
      continue
    else:
      chainfig = None
      chainfig = pppu.plot_posterior_chain(poslist, param, ifosNew, mcmcgr, \
                                           withhist=30)
      if chainfig:
        figname = output_fig(chainfig, puldir, param+fnamepostfix, ftypes)
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
  poststatstext = None
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
""" )

  if (swinj or hwinj) and parinj:
    poststatstext.append('<th>Inj. value</th>')

  poststatstext.append( \
"""\
    <th>IFO</th>
    <th>max. posterior</th>
    <th>mean</th>
    <th>median</th>
    <th>&sigma;</th>
    <th>%d%% confidence interval</th>
""" % int(opts.ci))

  if (swinj or hwinj) and parinj:
    poststatstext.append('<th>Inj. interval</th></tr>')
  else:
    poststatstext.append('</tr>')

  for j, param in enumerate(corrcoefsheader[0].split()):
    pl = param.lower()
    if 'post' in pl or 'logl' in pl or 'logw' in pl or 'logprior' in pl:
      continue

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

    if (swinj or hwinj) and parinj:
      try:
        poststatstext.append('<td rowspan="%d">%s</td>' % (len(ifosNew), \
 dispfunc(str(par[param.upper()]))))
      except:
        poststatstext.append('<td rowspan="%d">*</td>' % len(ifosNew))    

    for i, ifo in enumerate(ifosNew):
      poststatstext.append('<td class="%s">%s</td>' % (ifo, ifo) )

      poststatstext.append('<td>%s</td>' % dispfunc(str(maxLvals[i][param])))
      poststatstext.append('<td>%s</td>' % dispfunc(str(meanvals[i][param])))
      poststatstext.append('<td>%s</td>' % \
dispfunc(str(medvals[i][param])))
      poststatstext.append('<td>%s</td>' % \
dispfunc(str(stddevvals[i][param])))
      poststatstext.append('<td>(%s, %s)</td>' % \
(dispfunc(str(confidenceinterval[i][param][0])), dispfunc(str(confidenceinterval[i][param][1]))))

      if (swinj or hwinj) and parinj:
        poststatstext.append('<td>%.1f</td>' % confidenceregion[i][param])

      poststatstext.append('</tr>')

  poststatstext.append('</table>\n</div>')

  # output footer giving how the file was made
  pfootertext = []

  pfootertext.append( \
"""<div id="footer">
%s - %s <br>
<br>
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
  htmlptext.append(sdouttext)

  # output the posterior plots in an element
  htmlptext.append('<div class="wrapper">')
  htmlptext.append(cattext(posttabletext))
  htmlptext.append('</div>')

  # output the prior plots element
  if priortabletext:
    htmlptext.append('<div class="wrapper">')
    htmlptext.append(cattext(priortabletext))
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

  # close the shelf
  psrshelf.close()
