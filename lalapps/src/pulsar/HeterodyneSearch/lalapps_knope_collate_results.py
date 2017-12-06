# -*- coding: utf-8 -*-
#
#       lalapps_knope_collate_results.py
#
#       Copyright 2013, 2016
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

from __future__ import print_function, division

#standard library imports
import sys
import os
import re
import shutil
import datetime
import json
import argparse
from ConfigParser import ConfigParser
import subprocess as sp
import ast

#related third party imports
import numpy as np

import matplotlib
matplotlib.use("Agg")

#local application/library specific imports
from pylal import git_version
from lalapps.pulsarpputils import *
from lalapps.pulsarhtmlutils import *

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date


# create format for the output page
htmlpage = """
<!DOCTYPE html>
<html lang="en">

<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
  <meta name="description" content="Known pulsar search results"/>
  <meta charset="UTF-8">

  <title>Known pulsar search results</title>

  <link rel="stylesheet" type="text/css" href="{cssfile}"/>
</head>

<body>

<h1>Results table {nsources}</h1>

<!-- The table of pulsars and results -->
<div class="tablediv">
{resultstable}
</div>
<br />

<!-- Footnotes from the table -->
{footnotes}
<br />

<!-- a footer -->
<div id="footer">
{footer}
</div>

</body>
</html>
"""


# main function
if __name__=='__main__':
  description = """This script will collate results pages from multiple pulsars into a signal table."""
  epilog = """An example configuration file could contain the following:

# a section containing output information
[output]
path = path_to_output_directory  # the path to the output directory to contain the page [required]

# a section containing input directory information
[input]
path = path_to_input_directory   # the path to the directory in which the individual results page directories live [required]

# a section containing general information for the output table
[general]
detectors = ['H1', 'L1', 'Joint']   # a list of detctors (including 'Joint') whose results will be output
sort_value = name                   # the pulsar parameter on which to sort the results (default: 'name') [Allowed: 'name', 'freq', 'ra', 'dec', 'sdlim', 'hul', 'ell', 'sdrat', 'dist']
sort_direction = ascending          # sort in ascending or decending order (default: 'ascending')
results = ['h0ul', 'ell']            # a list of result parameters to output (default: ['h0ul'] - the h0 upper limit) [Allowed: 'h0ul', 'ell', 'sdrat', 'q22', 'bsn', 'bci', 'bcin', 'C21ul', 'C22ul']
parameters = ['f0rot', 'ra', 'dec'] # a list of pulsar parameters to output (default: ['f0rot'] - the pulsar's rotation frequency) [Allowed: 'f0rot', 'f0gw', 'f1rot', 'f1gw', 'sdlim', 'ra', 'dec', 'dist']
"""

  parser = argparse.ArgumentParser( description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter )
  parser.add_argument("inifile", help="The configuration (.ini) file")

  # parse input options
  opts = parser.parse_args()

  inifile = opts.inifile

  cp = ConfigParser()
  try:
    cp.read(inifile)
  except:
    print("Error... cannot parse configuration file '%s'" % inifile, file=sys.stderr)
    sys.exit(-1)

  # get the output path
  try:
    outpath = cp.get('output', 'path')
  except:
    print("Error... no output directory 'path' specified in [output] section.", file=sys.stderr)
    sys.exit(1)

  # try making directory
  if not os.access(outpath, os.W_OK) and not os.path.isdir(outpath): # check if directory already exists
    try:
      os.makedirs(outpath)
    except:
      print("Error... cannot make output directory '%s'." % outpath, file=sys.stderr)
      sys.exit(1)

  # get the directory containing all the individual pulsar result directories
  try:
    inpath = cp.get('input', 'path')
  except:
    print("Error... no input directory 'path' specified in [input] section.", file=sys.stderr)
    sys.exit(1)

  if not os.path.isdir(inpath):
    print("Error... input path '%s' does not exist." % inpath, file=sys.stderr)
    sys.exit(1)

  # check how to sort the results (e.g. by name)
  try:
    sorttype = cp.get('general', 'sort_value')
  except:
    sorttype = 'name' # default to sort by pulsar name

  # sorting can be by name (equivalently right ascension), frequency, declination, spin-down limit, h0 upper limit, ellipticity upper limit, spin-down ratio or distance
  if sorttype not in ['name', 'freq', 'ra', 'dec', 'sdlim', 'h0ul', 'ell', 'sdrat', 'dist']:
    print("Error... sorting can only be on either 'name', 'freq', 'ra', 'dec', 'sdlim', 'hul', 'ell', 'sdrat' or 'dist'.", file=sys.stderr)
    sys.exit(1)

  # get the direction of the sorting ('ascending' or 'descending')
  try:
    sortdirection = cp.get('general', 'sort_direction')
    if sortdirection == 'ascending':
      reverse = False
    elif sortdirection == 'descending':
      reverse = True
    else: # default to ascending
      reverse = False
  except:
    reverse = False # default to ascending

  # get list of detectors (and/or 'Joint') to use in the output table
  try:
    ifos = ast.literal_eval(cp.get('general', 'detectors'))
  except:
    print("Error... could not parse list of detectors to include.", file=sys.stderr)
    sys.exit(1)

  # get list of value to output (e.g. ['h0ul', 'ell'] for the h0 upper limit and ellipticity limit
  try:
    outputlims = [ol.upper() for ol in ast.literal_eval(cp.get('general', 'results'))] # covert to upper case
  except: # as default output the h0 upper limit
    outputlims = ['H0UL']

  # get a list of the pulsar parameters to output
  try:
    outputvals = [ov.upper() for ov in ast.literal_eval(cp.get('general', 'parameters'))] # convert to upper case
  except: # as default output the rotational frequency
    outputvals = ['F0ROT']

  # dictionary for output html and LaTeX pages
  htmlinput = {}
  latexinput = {}

  # get directories containing results pages
  resdirs = os.listdir(inpath)
  if len(resdirs) == 0:
    print("Error... the input directory '%s' contains no other directories." % inpath, file=sys.stderr)
    sys.exit(1)

  resultsdata = {} # dictionary to hold results data
  sourcedirs = [os.path.join(inpath, rd) for rd in resdirs if os.path.isdir(os.path.join(inpath, rd))]
  totalsources = len(sourcedirs)
  cursources = 0 # currently number of completed sources
  for d in sourcedirs:
    ld = os.listdir(d)
    jsonfile = None
    for fd in ld:
      if '.json' in fd:
        jsonfile = os.path.join(d, fd)
        break
    if jsonfile == None: # no file found, so move on to next directory
      continue

    # read in json file
    try:
      fp = open(jsonfile, 'r')
      pdata = json.load(fp)
      fp.close()
    except:
      print("Warning... could not read JSON file '%s'. Skipping this directory." % jsonfile)
      continue

    # check dictionary contains a 'Pulsar data' key
    if 'Pulsar data' not in pdata:
      print("Warning... no 'Pulsar data' field in JSON file '%s'. Skipping this directory." % jsonfile)
      continue

    # check dictionary contains 'PSR' key
    if 'PSR' not in pdata:
      print("Warning... no 'PSR' pulsar name field in JSON file '%s'. Skipping this directory." % jsonfile)
      continue

    # check that the request detectors are within the dictionary
    ifopresent = True
    for ifo in ifos:
      if ifo not in pdata:
        print("Warning... no key for detector '%s' in JSON file '%s'. Skipping this directory." % (ifo, jsonfile))
        ifopresent = False
        break
    if not ifopresent: continue

    # add data into dictionary
    resultsdata[pdata['PSR']] = pdata
    resultsdata[pdata['PSR']]['path'] = os.path.relpath(d, outpath) # relative path to result directory
    cursources += 1

  if len(resultsdata) == 0:
    print("Error... no reults were found!", file=sys.stderr)
    sys.exit(1)

  # perform sorting
  if 'freq' == sorttype: # sort by frequency
    sortlist = [(resultsdata[pname]['Pulsar data']['F0'], pname) for pname in resultsdata]
  elif 'dec' == sorttype: # sort by declination
    sortlist = [(resultsdata[pname]['Pulsar data']['DEC'], pname) for pname in resultsdata]
  elif 'dist' == sorttype: # sort by distance
    sortlist = [(resultsdata[pname]['Pulsar data']['DIST'], pname) for pname in resultsdata]
  elif 'sdlim' == sorttype: # sort by spin-down limit
    sortlist = [(resultsdata[pname]['Pulsar data']['spin-down limit'], pname) for pname in resultsdata]
  elif 'h0ul' == sorttype: # sort on h0 upper limit (get minimum h0 upper limit from all requested IFOs)
    sortlist = [(min([resultsdata[pname][ifo]['Upper limits']['H0'] for ifo in ifos]), pname) for pname in resultsdata]
  elif 'ell' == sorttype: # sort on ellipticity upper limit (get minimum upper limit from all requested IFOs)
    sortlist = [(min([resultsdata[pname][ifo]['Upper limits']['ELL'] for ifo in ifos]), pname) for pname in resultsdata]
  elif 'sdrat' == sorttype: # sort on ratio of result to spin-down limit (get minimum upper limit from all requested IFOs)
    sortlist = [(min([resultsdata[pname][ifo]['Upper limits']['spin-down ratio'] for ifo in ifos]), pname) for pname in resultsdata]
  else: # sort by name (equivalent to sorting by RA) by default
    sortlist = [(pname, pname) for pname in resultsdata]

  sortedlist = [p[1] for p in sorted(sortlist, reverse=reverse)] # return a list of sorted names

  nifos = len(ifos) # number of IFOs (including "Joint")
  numprepars = 1 + len(outputvals) # number of parameters in start of table
  numlims = len(outputlims) # number of limit parameter to output

  # check if outputing coherent vs incoherent Bayes factors (Bci or Bcin) - these will only be used for if 'Joint' is present in the detectors
  showbci = showbcin = 0
  if 'BCI' in outputlims:
    if 'Joint' in ifos:
      showbci = 1
    outputlims.pop(outputlims.index('BCI'))
    numlims -= 1
  if 'BCIN' in outputlims:
    if 'Joint' in ifos:
      showbcin = 1
    outputlims.pop(outputlims.index('BCIN'))
    numlims -= 1

  # start creating the results table
  restable = htmltable()
  ltable = latextable(ncolumns=(numprepars + numlims*nifos + showbci + showbcin))
  ltable.set_caption("Limits on the gravitational wave amplitude for known pulsars")

  # first row
  restable.addrow()
  ltable.addrow(underline=True)
  restable.adddata("", dataclass="bottomborder", header=True, colspan=numprepars)
  ltable.adddata("~", multicolumn=numprepars)
  for ifo in ifos:
    ncolspan = numlims
    if ifo == 'Joint':
      ncolspan += (showbci + showbcin)
    restable.adddata(ifo, dataclass=ifo, header=True, datastyle="text-align:center; border-left:1px solid #000; border-bottom:1px solid #000", colspan=ncolspan)
    ltable.adddata(ifo, multicolumn=ncolspan)

  # second row
  restable.addrow()
  ltable.addrow(underline=True)
  restable.adddata("Pulsar", dataclass="bottomborder", header=True)
  ltable.adddata("Pulsar")
  for prepar in outputvals:
    restable.adddata(paramhtmldict[prepar], dataclass="bottomborder", header=True)
    ltable.adddata(paramlatexdict[prepar])

  for ifo in ifos:
    dataclass = "leftborder" # class for first value
    datastyle = "border-bottom:1px solid #000" # style for first value
    for ol in outputlims:
      if ol[-2:] == 'UL': # if value is an upper limit add credible region
        cr = resultsdata.values()[0][ifo]['Upper limits']['credible region']
        restable.adddata(paramhtmldict[ol].format(cr), dataclass=dataclass, datastyle=datastyle, header=True)
        ltable.adddata(paramlatexdict[ol].format(cr))
      else:
        restable.adddata(paramhtmldict[ol], dataclass=dataclass, datastyle=datastyle, header=True)
        ltable.adddata(paramlatexdict[ol])
      dataclass = "bottomborder"
      datastyle = ""
    if ifo == 'Joint': # check whether to show additional coherent vs incoherent Bayes factor values
      if showbci:
        restable.adddata(paramhtmldict['BCI'], dataclass=dataclass, datastyle=datastyle, header=True)
        ltable.adddata(paramlatexdict['BCI'])
      if showbcin:
        restable.adddata(paramhtmldict['BCIN'], dataclass=dataclass, datastyle=datastyle, header=True)
        ltable.adddata(paramlatexdict['BCIN'])

  # a dictionary to convert between parameter names
  convdict = {'RA': 'RA', 'DEC': 'DEC', 'H0UL': 'H0', 'C21UL': 'C21', 'C22UL': 'C22', 'I21UL': 'I21', 'I31UL': 'I31', 'SDLIM': 'spin-down limit', 'SDRAT': 'spin-down ratio', 'BSN': 'Signal vs Noise', 'BCI': 'Coherent vs Incoherent', 'BCIN': 'Coherent vs Incoherent or Noise'}

  # checks for whegther footnotes are required
  dagger = False
  ddagger = False

  # loop through par files to produce plots
  for pname in sortedlist:
    restable.addrow()
    ltable.addrow()

    pulsar = resultsdata[pname] # the dictionary of information for this pulsar

    # output pulsar names and pulsar parameters (include link to results page)
    restable.adddata(atag(os.path.join(pulsar['path'], pname+'.html'), linktext=pname).text)
    ltable.adddata(re.sub('\-', '\\\\textminus', pname)) # substitute - ('hyphen') for \textminus in names

    for prepar in outputvals:
      htmlsdtag = ''
      latexsdtag = ''

      # set footnotes for pulsar's that have had a "corrected" spin-down
      if 'sdlim' in prepar:
        if pulsar['Pulsar data']['F1SD'] != None: # spin-down has been corrected for intrinsic motion effects
          htmlsdtag = htmltag('sup', tagtext="&dagger;").text
          latexsdtag = '$\dagger$'
          dagger = True
        elif pulsar['Pulsar data']['ASSOC'] != None:
          if 'GC' in pulsar['Pulsar data']['ASSOC']:
            htmlsdtag = htmltag('sup', tagtext="&Dagger;").text
            latexsdtag = '$\ddagger$'
            ddagger = True

      if prepar in convdict:
        pn = convdict[prepar]
      else:
        pn = prepar
      if pn in pulsar['Pulsar data']:
        preval = pulsar['Pulsar data'][pn]
        if preval == None:
          prevalhtml = '*'
          prevallatex = '*'
        else:
          disphtmlfunc = paramhtmldispfunc.__dict__[prepar]
          displatexfunc = paramlatexdispfunc.__dict__[prepar]
          prevalhtml = disphtmlfunc(preval)
          prevallatex = displatexfunc(preval)
      else:
        prevalhtml = '*'
        prevallatex = '*'

      restable.adddata(prevalhtml+htmlsdtag)
      ltable.adddata(prevallatex+latexsdtag)

    for ifo in ifos:
      dataclass = "leftborder"

      for limpar in outputlims:
        if limpar in convdict:
          ln = convdict[limpar]
        else:
          ln = limpar

        section = 'Upper limits' # generally values will be in the 'Upper limits' section
        if limpar == 'BSN': # if getting the Bayes factor look in 'Bayes factors' section
          section = 'Bayes factors'

        if ln in pulsar[ifo][section]:
          limval = pulsar[ifo][section][ln]

          if limval == None:
            limvalhtml = '*'
            linvallatex = '*'
          else:
            if limpar == 'BSN': # convert to log base 10
              limval = limval/np.log(10.)
            disphtmlfunc = paramhtmldispfunc.__dict__[limpar]
            displatexfunc = paramlatexdispfunc.__dict__[limpar]
            limvalhtml = disphtmlfunc(limval)
            limvallatex = displatexfunc(limval)
        else:
          limvalhtml = '*'
          limvallatex = '*'

        restable.adddata(limvalhtml, dataclass=dataclass)
        dataclass = ""
        ltable.adddata(limvallatex)

      # check if needing to output coherent vs incoherent Bayes factors
      if ifo == 'Joint':
        for bu, show in [('BCI', showbci), ('BCIN', showbcin)]:
          if show:
            bn = convdict[bu]
            if bn in pulsar[ifo]['Bayes factors']:
              bval = pulsar[ifo]['Bayes factors'][bn]
              if bval == None:
                bvalhtml = '*'
                bvallatex = '*'
              else:
                bval = bval/np.log(10.)
                disphtmlfunc = paramhtmldispfunc.__dict__[bu]
                displatexfunc = paramlatexdispfunc.__dict__[bu]
                bvalhtml = disphtmlfunc(bval)
                bvallatex = displatexfunc(bval)
            else:
              bvalhtml = '*'
              bvallatex = '*'

            restable.adddata(bvalhtml, dataclass=dataclass)
            dataclass = ""
            ltable.adddata(bvallatex)

  htmlinput['nsources'] = '(%d/%d sources)' % (cursources, totalsources)
  htmlinput['resultstable'] = restable.tabletext

  # add footnotes
  if dagger or ddagger:
    htmlfootnotes = htmltag('div', tagclass="footnotes")
    if dagger:
      htmlfootnotes += htmltag('sup', tagtext="&dagger;").text
      htmlfootnotes += " The spin-down limit was calculated using a spin-down "
      htmlfootnotes += atag("http://www.atnf.csiro.au/research/pulsar/psrcat/psrcat_help.html?type=normal&highlight=p1_i#p1_i", linktext="corrected").text
      htmlfootnotes += " for proper motion effects.\n"
      htmlfootnotes += "<br />\n"
      ltable.set_postamble("$\dagger$ The pulsar's spin-down is corrected for proper motion effects. \\\\\n")
    if ddagger:
      htmlfootnotes += htmltag('sup', tagtext="&ddagger;").text
      htmlfootnotes += " The spin-down limit was calculated using a characteristic spin-down age of 10<sup>9</sup> years.\n"
      ltable.set_postamble(ltable.postamble + "$\ddagger$ The pulsar's spin-down is calculated using a characteristic spin-down age of $10^9$ years.\n")
    htmlinput['footnotes'] = htmlfootnotes.text
  else:
    htmlinput['footnotes'] = ''

  # create CSS
  cssfile = os.path.join(outpath, 'table.css')
  fp = open(cssfile, 'w')
  fp.write(results_table_css)
  fp.close()

  htmlinput['cssfile'] = os.path.basename(cssfile)

  # get time/date for file creation
  now = datetime.datetime.now()

  # add footer containing author, date and command lines used for page
  htmlinput['footer'] = "{} - {}<br><br>Command lines used:<br>{}<br>{}<br>".format(__author__, now.strftime('%a %d %b %Y'), ' '.join(sys.argv), __version__)

  # create page
  try:
    htmlfile = os.path.join(outpath, 'index.html')
    fp = open(htmlfile, 'w')
    fp.write(htmlpage.format(**htmlinput))
    fp.close()
  except:
    print("Error... there was a problem outputting the html page.", file=sys.stderr)
    sys.exit(1)

  # output the LaTeX table
  try:
    latexfile = os.path.join(outpath, 'resultstable.tex')
    fp = open(latexfile, 'w')
    fp.write(ltable.tabletext)
    fp.close()
  except:
    print("Error... there was a problem outputting the LaTeX table.", file=sys.stderr)
    sys.exit(1)
