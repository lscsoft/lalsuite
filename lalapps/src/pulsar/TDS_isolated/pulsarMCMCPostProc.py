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

#related third party imports
import numpy as np

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

# convert a floating point number into a string in X.X x 10^Z format
def exp_str(f):
  s = '%.e' % f
  ssplit = s.split('e')
  return '%.1f&times;10<sup>%d</sup>' % (float(ssplit[0]), int(ssplit[1]))

# main function
if __name__=='__main__':
  from optparse import OptionParser
  
  description = \
"""This script is for creating a results output page and plots for the known
   pulsar analysis. It uses inputs from the old MCMC code
   lalapps_pulsar_parameter_estimation."""

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
     multiple MCMC directories from each IFO would be entered as follows:
     --ifo H1 --Bkfiles /home/user/finehetsH1
     --ifo L1 --Bkfiles /home/user/finehetsL1
     --mcmcdirs /home/user/pdfs1,/home/user/pdfs2,/home/user/pdfs3
  """
  parser.add_option("-m","--mcmcdirs", dest="mcmcdirs",\
                    help="A list "\
                    "of MCMC directories containing chain for all IFOs")
  
  # Turn on 2D kdes
  #parser.add_option("--twodkdeplots", action="store_true", default=False,
  #                  dest="twodkdeplots")
  
  # Turn on R convergence tests
  #parser.add_option("--RconvergenceTests", action="store_true", default=False,
  #                  dest="RconvergenceTests")
  
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
  
  # split MCMC directories
  mcmcdirs = (opts.mcmcdirs).split(',')
  
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
  
  # create header name style
  cssh1 = \
"""h1 {
  margin: 0 0 16px 0;
  padding: 0 0 16px 0;
  font-size: 20px;
  font-weight: bold;
  letter-spacing: 0px;
  border-bottom: 1px solid #999;
  font-family: Tahoma, Sans-Serif;
  text-shadow: 2px 2px 2px #ccc;
}
  
"""
  css.write(cssh1)
  
  # create par file style
  csspar = \
""".par {
  float: left;
  padding: 10px;
  border: 1px solid #999;
  box-shadow: 3px 3px 3px #888;
  font-family: monospace;
}
  
"""
  css.write(csspar)
 
  # create footer style
  cssfooter = \
"""#footer {
  border-top: 1px solid #999;
  padding: 15px;
  font-family: monospace;
}

"""
  css.write(cssfooter)
  
  # create link style
  csslink = \
""".par a:link{
  color: #000000
  text-decoration: none;
}

.par a:visited{
  color: #000000
  text-decoration: none;
}

"""
  css.write(csslink)
  
  csspostplot = \
""".posplot{
  height: 275px;
  padding: 8px;
  border: 1px solid #999;
  box-shadow: 2px 2px 2px #888;
}

"""
  css.write(csspostplot)
  
  # create a class for a rotated image
  cssrotate = \
""".rotated{
  padding: 8px;
  border: 1px solid #999;
  box-shadow: 2px 2px 2px #888;
  width: 275px;
  -moz-transform: rotate(90deg);
  -webkit-transform: rotate(90deg);
  -o-transform: rotate(90deg);
  -ms-transform: rotate(90deg);
}

"""
  css.write(cssrotate)

  css.close()
  
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
    
    # get time series and PSD plots
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
    
    # create html page
    htmlppage = os.path.join(puldir, 'index.html')
    try:
      htmlpout = open(htmlppage, 'w')
    except:
      print >> sys.stderr, 'Could not open html page %s' % htmlppage
     
    htmlpout.write('<!DOCTYPE html>\n')
    htmlpout.write('<html lang="en">\n\n') 
    
    htmlphead = \
"""
<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
  <meta name="description" content="PSR %s"/>
  <meta charset="UTF-8"> 

  <title>PSR %s</title>

  <link rel="stylesheet" type="text/css" href="../%s"/>
 
  </head>
""" % (pname, pname, cssname)
    htmlpout.write(htmlphead)
    
    htmlpout.write('<body>\n')
    
    # print out pulsar name to file
    htmlpout.write('<h1>PSR ' + pname + '</h1>\n')
    
    # copy par file to pulsar directory
    shutil.copy(parfile, os.path.join(puldir, pname + '.par'))
    
    # print out par file info
    htmlpout.write('<div class="par">\n')
    htmlpout.write('<a href="'+pname+'.par">\n')
    htmlpout.write('<table border="0">\n')
    for param in ['RAJ', 'DECJ', 'F0', 'F1', 'F2', 'PEPOCH']:
      pa = par[param]
      if pa:
        htmlpout.write('<tr>\n<td>'+param+'</td>\n')
        htmlpout.write('<td>'+str(pa)+'</td>\n</tr>')
    htmlpout.write('</table></a>\n')
    htmlpout.write('</div>\n\n')
    
    # get spin-down upper limit
    sdlim = pppu.spin_down_limit(f0, f1sd, dist)
    
    # get time series and PSD plots
    Bkdata = []
    plotpsds = False
    for i, ifo in enumerate(ifos):
      Bkdata.append(Bkfiles[i] + '/finehet_' + pname + '_' + ifo)
    Bkfigs, psdfigs = pppu.plot_Bks_ASDs( Bkdata, ifos, \
plotpsds=plotpsds, removeoutlier=8 )
    
    # output plots of time series and psd
    for i, ifo in enumerate(ifos):
      figname = 'Bk'+'_'+ifo+'.png'
      Bkplotpath = os.path.join(puldir, figname)
      Bkfigs[i].savefig(Bkplotpath)
      
      if plotpsds:
        figname = 'ASD'+'_'+ifo+'.png'
        ASDplotpath = os.path.join(puldir, figname)
        psdfigs[i].savefig(ASDplotpath)
    
    # loop over detectors
    ifosNew = ifos
    if nifos > 1:
      ifosNew.append('Joint')
    
    poslist = []
    
    for i, ifo in enumerate(ifosNew):
      # read in MCMC chains
      mcmc = []
      neffs = [] # number of effective samples for each chain
      for chaindir in mcmcdirs:
        cfile = chaindir + '/' + 'MCMCchain_' + pname + '_' + ifo
        
        if os.path.isfile(cfile):
          # load MCMC chain
          mcmcChain = np.loadtxt(cfile, comments='%')
          
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
        
        else:
          print >> sys.stderr, "File %s does not exist!" % cfile
          sys.exit(1)
        
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
      
      poslist.append(pos)
    
    # get posterior plots
    # h0
    bounds = [0, float("inf")]
    ul = 0.95
    histbins = 30
    h0Fig, ulvals = pppu.plot_posterior_hist( poslist, 'h0', ifosNew, \
                                              bounds, histbins, \
                                              ul, overplot=True )
    h0figname = 'h0post.png'
    h0plotpath = os.path.join(puldir, h0figname)
    h0Fig[0].savefig(h0plotpath)
    
    # convert h0 uls into ellipticity and spin-down ratio
    ell = []
    sdrat = []
    h0last = ulvals[-1]
    for h0ul in ulvals:
      ell.append(pppu.h0_to_ellipticity(h0ul, f0, dist))
      sdrat.append(h0ul/sdlim)
    
    sdlimstr = exp_str(sdlim)
    
    # print limits to html file
    htmlpout.write('<div>\n')
    htmlpout.write('<table border="1">\n')
    htmlpout.write('<tr>\n  <th>Detector</th>\n')
    htmlpout.write('  <th>h<sub>0</sub><sup>95%</sup></th>\n')
    htmlpout.write('  <th>&#949;</th>\n')
    htmlpout.write('  <th>spin-down ratio (%s)</th>\n</tr>\n' % sdlimstr)
    for i, ifo in enumerate(ifosNew):
      ulstr = exp_str(ulvals[i])
      ellstr = exp_str(ell[i])
      
      tableout = \
"""<tr>
  <td>%s</td>
  <td>%s</td>
  <td>%s</td>
  <td>%.2f</td>
</tr>
""" % (ifo, ulstr, ellstr, sdrat[i])
    htmlpout.write(tableout)
    htmlpout.write('</table>\n')
    htmlpout.write('</div>\n\n')
    
    # phi0
    bounds = [0, 2*math.pi]
    phi0Fig, ulvals = pppu.plot_posterior_hist( poslist, 'phi0', ifosNew, \
                                        bounds, histbins, \
                                        0, overplot=True )
    phi0figname = 'phi0post.png'
    phi0plotpath = os.path.join(puldir, phi0figname)
    phi0Fig[0].savefig(phi0plotpath)
    
    # cos(iota)
    bounds = [-1, 1]
    ciFig, ulvals = pppu.plot_posterior_hist( poslist, 'cosiota', ifosNew, \
                                      bounds, histbins, \
                                      0, overplot=True )
    cifigname = 'cipost.png'
    ciplotpath = os.path.join(puldir, cifigname)
    ciFig[0].savefig(ciplotpath)
    
    # psi
    bounds = [-math.pi/4, math.pi/4]
    psiFig, ulvals = pppu.plot_posterior_hist( poslist, 'psi', ifosNew, \
                                       bounds, histbins, \
                                       0, overplot=True )
    psifigname = 'psipost.png'
    psiplotpath = os.path.join(puldir, psifigname)
    psiFig[0].savefig(psiplotpath)
   
    # get h0 vs cos(iota) 2D posterior histrogram (if single detector of for
    # joint posterior)
    #bounds = [[0, 3*h0last], [-1, 1]]
    h0ciFig = pppu.plot_posterior_hist2D([poslist[-1]], ['h0', 'cosiota'], \
[ifosNew[-1]], bounds=None, nbins=[30, 30]) 
    h0cifigname = 'h0cipost.png'
    h0ciplotpath = os.path.join(puldir, h0cifigname)
    h0ciFig[0].savefig(h0ciplotpath)
    
    # produce output table of posterior plots
    htmlpout.write('<p>\n')
    htmlpout.write('<table border="0">\n')
    htmlpout.write('<tr>\n')
    htmlpout.write('<td><img class="posplot" src="%s"/></td>\n' % h0figname)
    htmlpout.write('<td></td>\n</tr>\n') # no entry
    htmlpout.write('<tr>\n')
    htmlpout.write('<td><img class="posplot" src="%s"/></td>\n' % h0cifigname)
    htmlpout.write('<td><img class="rotated" src="%s"/></td>\n' % cifigname)
    htmlpout.write('</tr>\n<tr>\n')
    htmlpout.write('<td><img class="posplot" src="%s"/></td>\n' % phi0figname)
    htmlpout.write('<td><img class="posplot" src="%s"/></td>\n' % psifigname)
    htmlpout.write('</table>\n')
    htmlpout.write('</p>\n\n')
    
    # output footer giving how the file was made
    htmlpout.write('<div id="footer">\n')
    htmlpout.write(__author__ + "<br>\n")
    htmlpout.write(__version__ + "<br>\n")
    htmlpout.write(__date__ + "<br>\n")
    #htmlpout.write(' '.join(args) + '\n')
    htmlpout.write('</div>\n\n')
    
    htmlpout.write('</body>\n\n')
    
    htmlpout.write('</html>\n')
    htmlpout.close()

