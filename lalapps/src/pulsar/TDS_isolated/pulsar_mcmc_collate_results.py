#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#       pulsar_mcmc_collate_results.py
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
from optparse import OptionParser
import subprocess as sp

#related third party imports
import numpy as np

import matplotlib
matplotlib.use("Agg")

#local application/library specific imports
from pylal import git_version
from lalapps import pulsarpputils as pppu

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

# create css file text
csstext = """
/* create body style */
body {
  font-family: Verdana, Geneva, "Trebuchet MS", sans-serif;
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


# convert a right ascension string in format 'hh:mm:ss.s' to a LaTeX string like
# H^h M^m S^s.ss
def ra_latexstr(ra):
  hms = ra.split(":")

  if len(hms) == 1:
    hms.append('0')
    hms.append('0')
  elif len(hms) == 2:
    hms.append('0')

  ss = ('%.2f' % float(hms[2])).split('.')

  return "$%s^{\\rm h}%s^{\\rm m}%s^{\\rm s}\!.%s$" % (hms[0].zfill(2), \
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


# convert a declination string in format 'dd:mm:ss.s' to a LaTeX string like
# dd^o mm' ss''.ss
def dec_latexstr(ra):
  dms = ra.split(":")

  if len(dms) == 1:
    dms.append('0')
    dms.append('0')
  elif len(dms) == 2:
    dms.append('0')

  ss = ('%.2f' % float(dms[2])).split('.')

  return "$%s^{\circ}%s'%s''\!.%s$" % ((re.sub('\+', '', dms[0])).zfill(2), \
dms[1].zfill(2), ss[0].zfill(2), ss[1].zfill(2))


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

    myfig.clf()
    matplotlib.pyplot.close(myfig)

  return fnameret


# html text to display for different parameter names
paramhtmltitledisp = {'RAJ': '&alpha;', 'DECJ': '&delta;', \
                 'RA': '&alpha;', 'DEC': '&delta;', \
                 'F0': 'f<sub>0</sub> (Hz)', 'F1': 'f<sub>1</sub> (Hz/s)', \
                 'F2': 'f<sub>2</sub> (Hz/s<sup>2</sup>)', \
                 'PEPOCH': 'epoch (MJD)', 'A1': 'a sin<it>i</i> (lt s)', \
                 'E': '<it>e</it>', 'EPS1': '&epsilon;<sub>1</sub>', \
                 'EPS2': '&epsilon;<sub>2</sub>', \
                 'T0': 'T<sub>0</sub> (MJD)', \
                 'TASC': 'T<sub>asc</sub> (MJD)', \
                 'OM': '&omega;<sub>0</sub> (deg)',
                 'PB': 'Period (days)', 'H0': 'h<sub>0</sub>', \
                 'COSIOTA': 'cos&iota;', 'PSI': '&psi; (rad)', \
                 'PHI0': '&phi;<sub>0</sub> (rad)', \
                 'PMRA': 'p.m. &alpha; (rad/s)', \
                 'PMDC': 'p.m. &delta; (rad/s)', \
                 'PMDEC': 'p.m. &delta; (rad/s)', \
                 'F0ROT': 'f<sub>rotation</sub> (Hz)', \
                 'F0GW': 'f<sub>GW</sub> (Hz)', \
                 'DIST': 'Distance (kpc)', \
                 'F1ROT': 'Spin-down<sub>rotation</sub> (Hz/s)', \
                 'F1GW': 'Spin-down<sub>GW</sub> (Hz/s)', \
                 'SDLIM': 'Spin-down limit', \
                 'ELL': '&#949;', \
                 'SDRAT': 'ratio', \
                 'H95': 'h<sub>0</sub><sup>95%</sup>', \
                 'H0PRIOR': 'h<sub>0</sub><sup>95%</sup> prior', \
                 'SDPOWRAT': 'power ratio (%)', \
                 'Q22': 'Q<sub>22</sub> (kg m<sup>2</sup>)'}

                 
# LaTeX text to display for different parameter names
paramlatextitledisp = {'RAJ': '$\\alpha$', 'RA': '$\\alpha$', \
                  'DECJ': '$\delta$', 'DEC': '$\delta$', \
                  'DIST': 'distance (kpc)', \
                  'ELL': '$\\varepsilon$', \
                  'H95': '$h_0^{95\%}$', \
                  'SDRAT': 'ratio', \
                  'SDLIM': 'spin-down limit', \
                  'F0ROT': '$\\nu_{\\rm rot}$ (Hz)', \
                  'F0GW': '$\\nu_{\\rm GW}$ (Hz)', \
                  'F1ROT': '$\dot{\\nu}_{\\rm rot}$ (Hz/s)', \
                  'F1GW': '$\dot{\\nu}_{\\rm GW}$ (Hz/s)', \
                  'SDPOWRAT': 'power ratio (\%)', \
                  'Q22': '$Q_{22}$ (kg\,m$^2$)', \
                  'H0PRIOR': '$h_0^{95\%}$ prior'}

                  
# a class containing function to output parameter vales in the appropriate
# format
class paramhtmlvaldisp:
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
  def ELL(f): return exp_str(float(f), 1)
  def H95(f): return exp_str(float(f), 1)
  def H0PRIOR(f): return exp_str(float(f), 1)
  def SDLIM(f): return exp_str(float(f), 1)
  def SDRAT(f): return '%.2f' % float(f)
  def DIST(f): return '%.1f' % float(f)
  def SDPOWRAT(f): return '%d' % int(f)
  def Q22(f): return exp_str(float(f), 1) # quadrupole moment
  def F0ROT(f): return '%.2f' % float(f)
  def F0GW(f): return '%.2f' % float(f)
  def F1ROT(f): return exp_str(float(f), 1)
  def F1GW(f): return exp_str(float(f), 1)
  def DEFAULTSTR(f): return f


# a class for outputting parameter values to the table
class paramlatexvaldisp:
  def ELL(f): return exp_latex_str(float(f), 1)
  def H95(f): return exp_latex_str(float(f), 1)
  def H0PRIOR(f): return exp_latex_str(float(f), 1)
  def SDLIM(f): return exp_latex_str(float(f), 1)
  def SDRAT(f): return '%.2f' % float(f)
  def RAJ(f): return ra_latexstr(f) # RA in string format
  def DECJ(f): return dec_latexstr(f) # dec in string format
  def RA(f): return ra_latexstr(f) # RA in string format
  def DEC(f): return dec_latexstr(f) # dec in string format
  def DIST(f): return '%.1f' % float(f)
  def SDPOWRAT(f): return '%d' % int(f)
  def Q22(f): return exp_latex_str(float(f), 1) # quadrupole moment
  def F0ROT(f): return '%.2f' % float(f)
  def F0GW(f): return '%.2f' % float(f)
  def F1ROT(f): return exp_latex_str(float(f), 1)
  def F1GW(f): return exp_latex_str(float(f), 1)
  def DEFAULTSTR(f): return f


# concatenate a list of strings (default to a new line between each)
def cattext(strlist, delim='\n'):
  return delim.join([strpart for strpart in strlist])


# main function
if __name__=='__main__':
  description = \
"""This script is for collating the results from multiple pulsars into a
table."""

  usage = "Usage: %prog [options]"

  parser = OptionParser( usage = usage, description = description,
                         version = __version__ )

  parser.add_option("-o", "--outpath", dest="outpath", help="The path for "
                    "the analysis output", metavar="DIR")

  parser.add_option("-z","--inpath", dest="inpath",\
                    help="The path to the directory containing the "
                    "individual pulsar directories pulsar", metavar="DIR")

  # get pulsar .par file directory
  parser.add_option("-p", "--parfile", dest="parfile", help="The "
                    "directory containing several, or an individual, "
                    "TEMPO-style pulsar parameter file used in the analysis. "
                    "[required]", metavar="PNAME.par", default=None)
                    
  # say if you want to compile the latex table (default to False)
  parser.add_option("-l", "--latex", dest="compilelatex", help="Set if you " \
                    "want to compile the latex results table",
                    action="store_true", default=False)

  # options for sorting the output (e.g. sort the by name [default], 
  # frequency, RA, dec)
  parser.add_option("-s", "--sort", dest="sorttype", help="Set the parameter " \
                    "that you want the output to be sorted by - one of either " \
                    "name (default), freq, ra or dec.", default="name") 
  
  # options for which IFOs to output (default is all IFO data used including the
  # joint analysis). "Joint" is assumed as an allowed IFO name. It the results
  # for a requested IFO don't exist it will be skipped.
  parser.add_option("-i", "--ifo", dest="ifolist", help="The individual IFOs "
                    "(or 'Joint' for results from a joint analysis) for which "
                    "you want results output (e.g. -i H1 -i L1 -i Joint). "
                    "[Default: use results from all available IFOs, including "
                    "any joint analysis]", 
                    action="append", default=None)
  
  # options for which calculated values to output (default is to output h0 upper
  # limit, ellipticity upper limit and spin-down ratio if available)
  parser.add_option("-u", "--outputlims", dest="outputlims", help="Say which upper "
                    "limit values that you want to output (can be 'h95', 'ell' "
                    "'sdrat', 'q22' or 'sdpow') e.g. -u h95 -u ell. [Default: output "
                    "h95, ell and sdrat]", action="append",
                    default=None)
  
  # options for which intrinsic values to output in the table (rotation frequency,
  # rotational fdot, distance and spin-down limit)
  parser.add_option("-n", "--outputvals", dest="outputvals", help="Say which pulsar "
                    "values you want to be output (can be 'f0rot', 'f0gw', 'f1rot', "
                    "'f1gw', 'dist', 'sdlim', 'ra', 'dec' or 'h0prior') e.g. -n "
                    "f0rot -n dist. [Default: output f0rot, f1rot, dist and sdlim]",
                    action="append", default=None)
  
  # output extra collated results plots
  # output histograms of the h95, ellipticity, spin-down ratio, and Q22 upper
  # limits (as in Figure 3 of http://adsabs.harvard.edu/abs/2010ApJ...713..671A)
  parser.add_option("-k", "--outputhists", dest="outputhists", help="Set to output "
                    "histograms of the upper limits.", action="store_true", 
                    default=False)
                    
  # output the h95 upper limits plotted against GW Frequency, and on top of the
  # estimates limit (calculated from the amplitude spectral densities)
  parser.add_option("-t", "--outputhulplot", dest="outputhulplot", help="Set to "
                    "output the h0 upper limits plot against frequency.", 
                    action="store_true", default=False)
                    
  # say whether to output any prior limits (e.g. S5 results on S6 plots) on the 
  # above two plots
  parser.add_option("-w", "--withprior", dest="withprior", help="Set to output any "
                    "prior results (if given) on the output histogram and/or upper "
                    "limit plot", action="store_true", default=False)
  
  # say if you want to output eps figure (default to False)
  parser.add_option("-e", "--epsout", dest="epsout", help="Set if wanting to " \
                    "output eps versions of figures", action="store_true",
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
  if not opts.__dict__['inpath']:
    print "Must specify in input directory containing the pulsars"
    parser.print_help()
    exit(-1)
  else:
    inpath = opts.inpath
    
  # create output directory if it doesn't exist
  if not os.path.isdir(outpath):
    try:
      os.mkdir(outpath)
    except:
      print >> sys.stderr, "Cannot create output directory %s" % outpath
      sys.exit(1)

  # if no particular outputs are specified default to h0 ul, ell ul and spin-down ratio ul
  if not opts.outputlims:
    outputlims = ['h0', 'ell', 'sdrat']
  else:
    outputlims = list(opts.outputlims)
  
  if not opts.outputvals:
    outputvals = ['f0rot', 'f1rot', 'dist', 'sdlim']
  else:
    outputvals = list(opts.outputvals)
  
  # check if parfile is a single file or a directory
  param_files = []
  if os.path.isfile(opts.parfile):
    param_files.append(opts.parfile)
  else:
    if os.path.isdir(opts.parfile):
      # list the parameter files in the directory
      pfiles = os.listdir(opts.parfile)

      for f in pfiles:
        if '.par' in f: # only use .par files
          param_files.append(opts.parfile + '/' + f)
    else:
      print >> sys.stderr, "No par file or directory specified!"
      sys.exit(1)

  sorttype = opts.sorttype
  
  ifolist = None
  if opts.ifolist:
    ifolist = list(opts.ifolist)
  
  # output plot values
  outputhulplot = opts.outputhulplot
  outputhists = opts.outputhists
  withprior = opts.withprior
  
  if opts.epsout:
    ftypes = ['png', 'eps']
  else:
    ftypes = ['png']
  
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

  # perform sorting
  names = []
  freqs = []
  ras = []
  decs = []
  for i, parfile in enumerate(param_files):
    try:
      par = pppu.psr_par(parfile)
    except:
      print >> sys.stderr, "Par file %s could not be opened!" % parfile
      # remove from list and move on to next pulsar
      param_files.remove(i)
      continue

    if par['PSRJ']:
      names.append(par['PSRJ'])
    else:
      names.append(None)
      
    if par['F0']:
      freqs.append(par['F0'])
    else:
      freqs.append(None)
      
    if par['RAJ']:
      ras.append(pppu.ra_to_rad(par['RAJ']))
    else:
      ras.append(None)
      
    if par['DECJ']:
      decs.append(pppu.dec_to_rad(par['DECJ']))
    else:
      decs.append(None)
  
  if "name" in sorttype:
    # index of names sorted in alphanumeric order
    sortidx = sorted(range(len(names)), key=names.__getitem__)
  elif "freq" in sorttype:
    # freqs in ascending order
    sortidx = sorted(range(len(freqs)), key=freqs.__getitem__)
  elif "ra" in sorttype:
    # RAs in ascending order
    sortidx = sorted(range(len(ras)), key=ras.__getitem__)
  elif "dec" in sorttype:
    # decs in ascending order
    sortidx = sorted(range(len(decs)), key=decs.__getitem__)
  
  if outputhists or outputhulplot:
    h95 = {}
    ell = {}
    q22 = {}
    sdrat = {}
    ulestbot = {}
    ulesttop = {}
    f0gw = []
    hifos = []
  
  if withprior:
    h0prior = []
    ellprior = []
    q22prior = []
    sdratprior = []
    f0gwprior = []
  
  firstpsr = True
  
  # loop through par files to produce plots
  for idx in sortidx:
    print 'Results for pulsar ' + names[idx]

    skip = False
    
    # create directory containing pulsar data (assume it has been created 
    # using in the PSRJ name (as will be the case with pulsar_mcmc_create_results_page.py
    puldir = os.path.join(inpath, names[idx])
    if not os.path.isdir(puldir):
      print >> sys.stderr, "No directory for PSR%s, skip this pulsar" % names[idx]
      continue

    # get the shelved database of information 
    try:
      psrshelf = shelve.open(os.path.join(puldir, names[idx]+'.db'))
    except:
      print >> sys.stderr, "No database information for PSR%s. Skip this pulsar" % names[idx]
      psrshelf.close()
      continue
    
    sifos = sintrinsicsd = sagebasedsd = None
      
    # get some of the shelved data
    try:
      sifos = psrshelf['ifos'] # get IFOs used
      sintrinsicsd = psrshelf['intrinsicsd'] # get whether fdot has been corrected for proper motion
      sagebasedsd = psrshelf['agebasedsd'] # get whether fdot has been calculated used a characteristic age of 10^9 years
    except:
      continue
    
    # see which IFOs are being used
    if not ifolist:
      ifolist = sifos
    else:
      ifolisttmp = []
      for ifo in ifolist:
        if ifo in sifos:
          ifolisttmp.append(ifo)
          
      if not ifolisttmp:
        print >> sys.stderr, "No data for the specified IFOs: %s!" % ', '.join(ifolist)
        sys.exit(1)
      else:
        ifolist = ifolisttmp
    
    nifos = len(ifolist) # number of IFOs (including "Joint")
    numprepars = 1 + len(outputvals) # number of parameters in start of table    
    numlims = len(outputlims) # number of limit parameter to output
    
    # if on first available pulsars create the output page table header
    if firstpsr:
      htmltext.append('<tr>\n<th class="bottomborder" colspan="%d"></th>' % \
numprepars)
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
\\tablehead{\multicolumn{%d}{c}{~} \
""" % (' '.join(['c']*(numprepars + numlims*nifos)), numprepars) )

      for ifo in ifolist:
        htmltext.append( \
"""<th class="%s" colspan="%d" style="text-align:center; border-left:1px \
solid #000; border-bottom:1px solid #000">%s</th>""" % (ifo, numlims, ifo))
        latextext.append(' & \multicolumn{%d}{c}{%s}' % (numlims, ifo))

      htmltext.append('</tr><tr><th class="bottomborder">Pulsar</th>')
      latextext.append("\\\\ \colhead{Pulsar} ")
      
      for prepar in outputvals:
        htmltext.append( '<th class="bottomborder">%s</th>' % paramhtmltitledisp[prepar.upper()] )
        latextext.append( "& \colhead{%s} " % paramlatextitledisp[prepar.upper()] )

      for ifo in ifolist:
        cl = '<th class="leftborder" style="border-bottom:1px solid #000">'
        for limpar in outputlims:
          htmltext.append("%s%s</th>" % (cl, paramhtmltitledisp[limpar.upper()]))
          cl = '<th class="bottomborder">'
          
          latextext.append("& \colhead{%s}" % paramlatextitledisp[limpar.upper()])

      htmltext.append('</tr>')
      latextext.append('}\n\startdata')

      if outputhists or outputhulplot:
        hifos = list(ifolist)
      
      firstpsr = False
      
    else: # check that the number of IFOS
      if nifoprev != nifos:
        print >> sys.stderr, "Number of IFOs for pulsar %s is not the same as previous pulsar. Skip." % names[idx]
        psrshelf.close()
        continue
      else: # check the ifos are the same
        for ifo in ifolist:
          if ifo not in ifosprev:
            print >> sys.stderr, "An IFOs for pulsar %s is not the same as for the previous pulsar. Skip." % names[idx]
            psrshelf.close()
            continue
    
    ifosprev = list(ifolist)
    nifoprev = nifos
    
    # print out info to the main table
    htmlptext = []
    latexptext = []
    
    # output pulsar names and pulsar parameters (include link to results page)
    htmlptext.append('<tr><td><a href="%s/">%s</a></td>' % (names[idx], names[idx]))
    latexptext.append("%s " % re.sub('\-', '\\\\textminus', names[idx]))
    for prepar in outputvals:
      htmlsdtag = ''
      latexsdtag = ''
      if 'sdlim' in prepar:
        if sintrinsicsd:
          htmlsdtag = '<sup>&dagger;</sup>'
          latexsdtag = '$^{\dagger}$'
        elif sagebasedsd:
          htmlsdtag = '<sup>&Dagger;</sup>'
          latexsdtag = '$^{\ddagger}$'
      
      try:
        preval = psrshelf[prepar]
      except:
        print >> sys.stderr, "Error getting %s for PSR%s. Skipping pulsar." % (prepar, names[idx])        
        skip = True
        break
      
      if not preval:
        disphtmlfunc = paramhtmlvaldisp.__dict__['DEFAULTSTR']
        displatexfunc = paramlatexvaldisp.__dict__['DEFAULTSTR']
        preval = '*'
      else:
        disphtmlfunc = paramhtmlvaldisp.__dict__[prepar.upper()]
        displatexfunc = paramlatexvaldisp.__dict__[prepar.upper()]
      
      if skip:
        break
      
      htmlptext.append("<td>%s%s</td>" % (disphtmlfunc(str(preval)), htmlsdtag))
      latexptext.append(" & %s%s" % (displatexfunc(str(preval)), latexsdtag))
    
    for i, ifo in enumerate(ifolist):
      cl = '<td class="leftborder">'
      
      for limpar in outputlims:
        try:
          limvals = psrshelf[limpar]
        except:
          print >> sys.stderr, "Error getting %s for PSR%s. Skipping pulsar." % (limpar, names[idx])        
          skip = True
          break
        
        limval = limvals[ifo]
        if not limval:
          disphtmlfunc = paramhtmlvaldisp.__dict__['DEFAULTSTR']
          displatexfunc = paramlatexvaldisp.__dict__['DEFAULTSTR']
          limval = '*'
        else:
          disphtmlfunc = paramhtmlvaldisp.__dict__[limpar.upper()]
          displatexfunc = paramlatexvaldisp.__dict__[limpar.upper()]
        
        htmlptext.append('%s%s</td>' % (cl, disphtmlfunc(str(limval))))
        latexptext.append('& %s ' % displatexfunc(str(limval)))
        
        cl = '<td>'
        
      if skip:
        break
    
    if not skip:
      # get values to plot
      if outputhulplot or outputhists:
        # get upper an lower amplitude upper limit estimates
        try:
          ulb = psrshelf['ulestimatebot']
          ult = psrshelf['ulestimatetop']
          h95s = psrshelf['h95']
          ells = psrshelf['ell']
          sdrats = psrshelf['sdrat']
          q22s = psrshelf['q22']
          f0 = psrshelf['f0rot']
        except:
          ulb = ult = h95s = f0 = ells = sdrats = q22s = None
 
        if ulb and ult and h95s and f0 and ells and q22s and sdrats:
          f0gw.append(2.*f0)

          # create dictionary of lists
          if not h95:
            for ifo in hifos:
              h95[ifo] = [h95s[ifo]]
              ell[ifo] = [ells[ifo]]
              q22[ifo] = [q22s[ifo]]
              sdrat[ifo] = [sdrats[ifo]]
              ulestbot[ifo] = [ulb[ifo]]
              ulesttop[ifo] = [ult[ifo]]
          else:
            for ifo in hifos:
              h95[ifo].append(h95s[ifo])
              ell[ifo].append(ells[ifo])
              q22[ifo].append(q22s[ifo])
              sdrat[ifo].append(sdrats[ifo])
              ulestbot[ifo].append(ulb[ifo])
              ulesttop[ifo].append(ult[ifo])
      
      # get prior upper limit on h0
      if withprior:
        try:
          h0p = psrshelf['h0prior']
          dist = psrshelf['dist']
          f0 = psrshelf['f0rot']
          sdlim = psrshelf['sdlim']
        except:
          h0p = dist = f0 = sdlim = None
        
        if h0p and dist and f0 and sdlim:
          h0prior.append(h0p)
          ellprior.append(pppu.h0_to_ellipticity(h0p, f0, dist))
          q22prior.append(pppu.h0_to_quadrupole(h0p, f0, dist))
          f0gwprior.append(2.*f0)
          sdratprior.append(h0p/sdlim)
      
      # add pulsar to table
      htmltext.append(cattext(htmlptext))
      latextext.append(cattext(latexptext))

    psrshelf.close()
 
  
  htmltext.append('</table>')

  htmltext.append( \
"""
<p>
  <sup>&dagger;</sup> The spin-down limit was calculated using a spin-down
  <a \
href="http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=\
normal&highlight=p1_i#p1_i">corrected</a> for proper motion effects.
<br>
  <sup>&Dagger;</sup> The spin-down limit was calculated using a characteristic spin-down
  age of 10<sup>9</sup> years.
</p>
""" )

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
  latextext.append("\\tablenotetext{\dagger}{The pulsar's spin-down is \
corrected for proper motion effects.}")
  latextext.append("\\tablenotetext{\ddagger}{The pulsar's spin-down is \
calculated using a characteristic spin-down ago of $10^9$ years.}")
  latextext.append('\\end{deluxetable}')
  latextext.append('\clearpage')
  latextext.append('\\end{landscape}')
  latextext.append('\\end{document}')
  latexout.write(cattext(latextext))
  latexout.close()

  # compile the LaTeX table
  if opts.compilelatex:
    sp.call(['latex', latexpage])
    sp.call(['latex', latexpage])
    sp.call(['dvips', re.sub('tex', 'dvi', latexpage)])
    sp.call(['ps2pdf', re.sub('tex', 'ps', latexpage)])

  # create output histograms
  if outputhists:
    prevlims = None
    if withprior:
      prevlims = h0prior
   
    h95histfigs = pppu.plot_limits_hist(h95, 'h95', hifos, prevlims=prevlims)
    
    prevlims = None
    if withprior:
      prevlims = ellprior
    
    ellhistfigs = pppu.plot_limits_hist(ell, 'ell', hifos, prevlims=prevlims)
      
    prevlims = None
    if withprior:
      prevlims = sdratprior
    
    sdrathistfigs = pppu.plot_limits_hist(sdrat, 'sdrat', hifos, prevlims=prevlims)
      
    prevlims = None
    if withprior:
      prevlims = q22prior
    
    q22histfigs = pppu.plot_limits_hist(q22, 'q22', hifos, prevlims=prevlims)
    
    for i, ifo in enumerate(hifos):
      h95figname = output_fig(h95histfigs[i], outpath, 'h95hist_'+ifo, ftypes)
      ellfigname = output_fig(ellhistfigs[i], outpath, 'ellhist_'+ifo, ftypes)
      sdratfigname = output_fig(sdrathistfigs[i], outpath, 'sdrathist_'+ifo, ftypes)
      q22figname = output_fig(q22histfigs[i], outpath, 'q22hist_'+ifo, ftypes)
      
  # create output upper limit plot
  if outputhulplot:
    prevlims = prevf0 = None
    if withprior:
      prevlims = h0prior
      prevf0 = f0gwprior
    
    freqrange = [10, 1500]
    hulfigs = pppu.plot_h0_lims(h95, f0gw, hifos, freqrange, ulesttop, ulestbot, prevlims, prevf0)
    
    # output the plots
    for i, ifo in enumerate(hifos):
      hulfigname = output_fig(hulfigs[i], outpath, 'hul_'+ifo, ftypes)
      
