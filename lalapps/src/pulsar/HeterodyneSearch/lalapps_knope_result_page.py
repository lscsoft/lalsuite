# -*- coding: utf-8 -*-
#
#       lalapps_knope_result_page.py
#
#       Copyright 2015
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

# Script from creating results pages for pulsars from the known pulsar search

from __future__ import print_function, division

import argparse
from ConfigParser import ConfigParser
import sys
import ast
import numpy as np
import re
import urllib2
import copy
import os
import fnmatch
import ast
import datetime
import json
from scipy import stats
import h5py

import matplotlib
matplotlib.use("Agg")

from lalapps.pulsarpputils import *
from lalapps.pulsarhtmlutils import *
from pylal import bayespputils as bppu
from pylal import git_version

__author__="Matthew Pitkin <matthew.pitkin@ligo.org>"
__version__= "git id %s"%git_version.id
__date__= git_version.date

# try importing scotchcorner
try:
  from scotchcorner import scotchcorner
except ImportError:
  print("Could no import scotchcorner: make sure scotchcorner is installed (e.g. 'pip install scotchcorner') and in the PYTHONPATH", file=sys.stderr)
  sys.exit(1)

# create format for the output page
htmlpage = """
<!DOCTYPE html>
<html lang="en">

<head>
  <meta http-equiv="Content-Type" content="text/html; charset=utf-8"/>
  <meta name="description" content="PSR {psrname}"/>
  <meta charset="UTF-8">

  <title>{title}</title>

  <link rel="stylesheet" type="text/css" href="{cssfile}"/>

</head>

<body>

<!- add in javascript to allow text toggling -->
<script language="javascript">
function toggle(id) {{
  var ele = document.getElementById(id);
  if(ele.style.display == "block"){{
    ele.style.display = "none";
  }}
  else
    ele.style.display = "block";
}}
</script>

<!-- Page title -->
<h1>{h1title}</h1>

<!-- Links to the parts of the page -->
<div class="pagelinks">
{linkstable}
</div>

<!-- table of pulsar parameters from EM data -->
<table>
<tr>
<td>
<div class="pulsartable" id="pulsartable">
{pulsartable}
</div>
</td>

<td>
<!-- table of derived gravitational wave upper limits, SNRs and evidence ratios -->
<div class="limitstable" id="limitstable">
{limitstable}
</div>
</td>
</tr>
</table>

<br />

<!-- plots of 1D and 2D marginalised posteriors for GW amplitude and phase parameters -->
<div class="selectedposteriors" id="selectedposteriors">
{selectedposteriors}
</div>

<br />

<!-- corner plot of all parameters (hidden by default) -->
<div class="jointposteriors" id="jointposteriors">
{jointposteriors}
</div>

<br />

<!-- (if background runs have been produced): if running a multi-detector analysis the plot
     coherent vs. incoherent Bayes factors (with SNR giving colour of points) of background
     and foreground, or if single detector plot signal vs. noise Bayes factor against SNR. -->
<div class="evidenceplots" id="evidenceplots">
{evidenceplots}
</div>

<br />

<!-- plots of pre-processed data (heterodyned/spectrally interpolated time series) -->
<div class="dataplots" id="dataplots">
{dataplots}
</div>

<br />

<!-- plots of the posterior samples for all parameters -->
<div class="posteriorsamples" id="posteriorsamples">
{posteriorsamples}
</div>

<br />

<!-- table of posterior statistics: means, max. likelihoods, confidence intervals -->
<div class="posteriorstats" id="posteriorstats">
{posteriorstats}
</div>

<br />

<!-- a footer -->
<div id="footer">
{footer}
</div>

</body>
</html>
"""


def get_atnf_info(psr):
  """
  Get the pulsar (psr) distance (DIST in kpc), proper motion corrected age (AGE_I) and any association
  (ASSOC e.g. GC) from the ATNF catalogue.
  """

  psrname = re.sub('\+', '%2B', psr) # switch '+' for unicode character

  atnfversion = '1.54' # the latest ATNF version
  atnfurl = 'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=' + atnfversion
  atnfurl += '&Dist=Dist&Assoc=Assoc&Age_i=Age_i' # set parameters to get
  atnfurl += '&startUserDefined=true&c1_val=&c2_val=&c3_val=&c4_val=&sort_attr=jname&sort_order=asc&condition=&pulsar_names=' + psrname
  atnfurl += '&ephemeris=selected&submit_ephemeris=Get+Ephemeris&coords_unit=raj%2Fdecj&radius=&coords_1=&coords_2='
  atnfurl += '&style=Long+with+last+digit+error&no_value=*&fsize=3&x_axis=&x_scale=linear&y_axis=&y_scale=linear&state=query'

  try:
    urldat = urllib2.urlopen(atnfurl).read() # read ATNF url
    predat = re.search(r'<pre[^>]*>([^<]+)</pre>', urldat) # to extract data within pre environment (without using BeautifulSoup) see e.g. http://stackoverflow.com/a/3369000/1862861 and http://stackoverflow.com/a/20046030/1862861
    pdat = predat.group(1).strip().split('\n') # remove preceeding and trailing new lines and split lines
  except:
    print("Warning... could not get information from ATNF pulsar catalogue.")
    return None

  # check whether information could be found and get distance, age and association from data
  dist = None
  age = None
  assoc = None
  for line in pdat:
    if 'WARNING' in line or 'not in catalogue' in line:
      return None
    vals = line.split()
    if 'DIST' in vals[0]:
      dist = float(vals[1])
    if 'AGE_I' in vals[0]:
      age = float(vals[1])
    if 'ASSOC' in vals[0]:
      assoc = vals[1]

  return (dist, age, assoc)


def set_spin_down(age, assoc, f0, f1):
  """
  Set the spin-down of the source based on the instrinsic age (age) corrected for any proper motion/
  globular cluster accelation if available, or if not give AND the pulsar is in a globular cluster base the
  spin-down on assuming an age of 10^9 years. Otherwise just return the unadjusted spin-down.
  """

  if age != None and age > 0.:
    return -f0/(2.* age * 365.25 * 86400.)
  elif assoc != None:
    if 'GC' in assoc: # check if a globular cluster pulsar
      return -f0/(2. * 1.e9 * 365.25 * 86400.)
    else:
      return f1
  else:
    return f1


def create_psr_table(par):
  """
  Create a html table of some information from the pulsar parameter file.
  """
  table = htmltable() # create table class

  # list of parameters to display (in this order)
  paramdisplist = ['RAJ', 'DECJ', 'F0', 'F1', 'F2', 'PEPOCH', 'A1', 'E', 'EPS1', 'EPS2', 'OM', 'T0', 'TASC', 'PB']

  for param in paramdisplist:
    pa = par[param]
    if pa != None:
      table.addrow()
      dispfunc = paramhtmldispfunc.__dict__[param]
      table.adddata(paramhtmldict[param])
      table.adddata(dispfunc(str(pa)))

  return table.tabletext


class create_data_table:
  def __init__(self, datafiles, outputdir, figformats=['png'], asdtime=86400, harmonics=[2.]):
    """
    Initialise with a dictionary keyed in detector names containing paths to the equivalent pre-processed
    data file(s) for that detector (if harmonics contains more than one frequency harmonic this there should
    be a list of files, with one for each harmonic and given in the same order as the harmonics lists).
    figformats gives a list of the output figure formats (defaulting to .png), and asd time is the time to use
    for each FFT when creating the spectrograms/ASDs, which defaults to 86400 seconds (1 day).
    """

    self._datatable = htmltag('h2', 'Analysis statistics', newline=True).text # the output html table

    if not isinstance(harmonics, list):
      self._harmonics = list(harmonics)
    else:
      self._harmonics = harmonics
    self._nharmonics = len(self._harmonics)
    self._ifos = list(datafiles.keys()) # get list of detectors
    self._datafiles = {} # dictionary (keyed to harmonics) of dictionaries (keyed to detector) of data files

    for i, ifo in enumerate(datafiles):
      if isinstance(datafiles[ifo], list):
        if len(datafiles[ifo]) != self._nharmonics:
          print("Error... number of pre-processed data files is not the same as the number of harmonics.", file=sys.stderr)
          sys.exit(1)
        else:
          if i == 0:
            for h in self._harmonics: self._datafiles[h] = {}

          for j, df in enumerate(datafiles[ifo]):
            if not os.path.isfile(df):
              print("Error... pre-processed data file '%s' for '%s' does not exist." % (df, ifo), file=sys.stderr)
              sys.exit(1)
            else:
              self._datafiles[self._harmonics[j]][ifo] = df
      else:
        if self._nharmonics == 1:
          self._datafiles[self._harmonics[0]] = datafiles
          for df in list(self._datafiles[self._harmonics[0]].values()):
            if not os.path.isfile(df):
              print("Error... pre-processed data file '%s' for '%s' does not exist." % (df, ifo), file=sys.stderr)
              sys.exit(1)
          break
        else:
          print("Error... number of pre-processed data files is not the same as the number of harmonics.", file=sys.stderr)
          sys.exit(1)

    self._outputdir = outputdir
    if not os.path.isdir(self._outputdir):
      print("Error... output path '%s' for data plots does not exist" % self._outputdir, file=sys.stderr)
      sys.exit(1)

    self._figformats = figformats # figure output formats
    if 'png' not in figformats:
      print("Error... must include 'png' in the output figure formats", file=sys.stderr)
      sys.exit(1)

    # get plots of absolute values, amplitude spectral density and spectrograms of the pre-processed data
    self._asdtime = asdtime
    self._asds = {}
    self.create_plots()

    self._starts = {} # start times of data files
    self._ends = {}   # end times of data files
    self._length = {} # length (seconds) of data files
    self._duty_factor = {} # duty factors of data

    # get analysis data statistics (start time, end time, duty factor)
    self.get_stats()

    # create table of results for all detectors
    self.create_table()

  def create_plots(self):
    # create the plots of the absolute values, amplitude spectral density and spectrograms of the
    # pre-processed data.
    self._Bkplots = {}
    self._asdplots = {}
    self._fscanplots = {}

    for h in self._harmonics:
      self._Bkplots[h] = {}
      self._asdplots[h] = {}
      self._fscanplots[h] = {}
      self._asds[h] = {}

      Bkfigs, asdfigs, fscanfigs, asdlist, _ = plot_Bks_ASDs(self._datafiles[h], delt=self._asdtime, plotpsds=True, plotfscan=True, removeoutlier=5.)

      if self._nharmonics == 1:
        harmtext = '' # don't add any harmonics suffix to the file name
      else:
        harmtext = '_%d' % int(h)

      # output the plots
      for i, ifo in enumerate(self._ifos):
        self._asds[h][ifo] = asdlist[i] # set the estimate of the amplitude spectral density for that IFO

        for ftype in self._figformats:
          for fig, fignameprefix in zip([Bkfigs[i], fscanfigs[i], asdfigs[i]], ['Bk', 'fscan', 'ASD']):
            # try outputting figures
            figname = os.path.join(outdir, '%s_%s%s.%s' % (fignameprefix, ifo, harmtext, ftype))
            try:
              fig.savefig(figname)
            except:
              print("Error... problem creating figure '%s'." % figname, file=sys.stderr)
              sys.exit(1)

          if ftype == 'png': # add figure files
            self._Bkplots[h][ifo] = os.path.join(outdir, 'Bk_%s%s.%s' % (ifo, harmtext, ftype))
            self._asdplots[h][ifo] = os.path.join(outdir, 'ASD_%s%s.%s' % (ifo, harmtext, ftype))
            self._fscanplots[h][ifo] = os.path.join(outdir, 'fscan_%s%s.%s' % (ifo, harmtext, ftype))

  @property
  def Bkplots(self):
    # return dictionary of plots of absolute values of the pre-processed data for each detector
    return self._Bkplots

  @property
  def asdplots(self):
    # return dictionary of plot of the amplitude specatral densities for each detector
    return self._asdplots

  @property
  def fscanplots(self):
    # return dictionary of spectrogram plots of the pre-processed data for each detector
    return self._fscanplots

  @property
  def asds(self):
    # return dictionary of noise amplitude spectral density estimates for each detector
    return self._asds

  def get_stats(self):
    # get analysis statistics (data start, end and duty cycle) for each detector
    for i, ifo in enumerate(self._ifos):
      # use np.loadtxt to open file and get start and end times
      try:
        bkd = np.loadtxt(self._datafiles[self._harmonics[0]][ifo], comments=['#', '%'])
      except:
        print("Error... could not load file '%s'." % self._datafiles[self._harmonics[0]][ifo], file=sys.stderr)
        sys.exit(1)

      self._starts[ifo] = bkd[0,0] # start time
      self._ends[ifo] = bkd[-1,0]  # end time
      dt = np.min(np.diff(bkd[:,0])) # minimum time step in data
      self._length[ifo] = float(len(bkd))*dt
      self._duty_factor[ifo] = 100.*self._length[ifo]/(self._ends[ifo]-self._starts[ifo])

  def analysis_stats_td(self, ifo):
    # return the inside of the html <td> element giving the statistics for a particular detector
    table = htmltable()
    table.addrow()
    table.adddata('&nbsp;')
    table.adddata(ifo, dataclass=ifo)
    tdvals = [('Start (GPS)', str(int(self._starts[ifo]))), ('End (GPS)', str(int(self._ends[ifo]))), ('Length (sec)', str(int(self._length[ifo]))), ('Duty factor (%)', '%.lf' % self._duty_factor[ifo])]
    for item in tdvals:
      table.addrow()
      table.adddata(item[0])
      table.adddata(item[1])

    return table.tabletext

  def create_table(self):
    table = htmltable()

    for ifo in self._ifos:
      for h in self._harmonics:
        table.addrow()
        if self._nharmonics != 1:
          table.adddata('Data for %df' % int(h), colspan=2, header=True)
          table.addrow()

        table.adddata(self.analysis_stats_td(ifo), datastyle="text-align: center;")
        plotlink = atag(os.path.basename(self._Bkplots[h][ifo]), '<img class="dataplot" src="{}"/>'.format(os.path.basename(self._Bkplots[h][ifo]))).text
        table.adddata(plotlink)
        table.addrow()
        plotlink = atag(os.path.basename(self._asdplots[h][ifo]), '<img class="asdplot" src="{}"/>'.format(os.path.basename(self._asdplots[h][ifo]))).text
        table.adddata(plotlink)
        plotlink = atag(os.path.basename(self._fscanplots[h][ifo]), '<img class="dataplot" src="{}"/>'.format(os.path.basename(self._fscanplots[h][ifo]))).text
        table.adddata(plotlink)

    self._datatable += table.tabletext

  def __str__(self): # string method
    return self._datatable


class posteriors:
  """
  Get sample posteriors and created a set of functions for outputting tables, plots and posterior statistics
  """
  def __init__(self, postfiles, outputdir, harmonics=[2], modeltype='waveform', biaxial=False, usegwphase=False, injection=None):
    """
    Initialise with a dictionary keyed in detector names containing paths to the equivalent posterior samples
    file for that detector.
    """
    self._outputdir = outputdir
    if not os.path.isdir(self._outputdir):
      print("Error... output path '%s' for data plots does not exist" % self._outputdir, file=sys.stderr)
      sys.exit(1)

    self._ifos = list(postfiles.keys()) # get list of detectors
    self._postfiles = postfiles
    self._posteriors = {}             # dictionary of posterior objects
    self._posterior_stats = {}        # dictionary if posteriors statistics
    self._signal_evidence = {}        # dictionary of signal evidence values
    self._noise_evidence = {}         # dictionary of noise evidence values
    self._maxL = {}                   # dictionary of maximum log likelihood values
    self._Bsn = {}                    # dictionary of signal vs noise Bayes factors
    self._signal_evidence = {}        # dictionary of signal evidences
    self._noise_evidence = {}         # dictionary of noise evidences
    self._Bci = None                  # coherent versus incoherent Bayes factor
    self._Bcin = None                 # coherent versus incoherent or noise Bayes factor
    self._optimal_snrs = {}           # dictionary of optimal matched filter SNRs
    self._parameters = []             # list of the source parameters in the posterior files
    self._injection = injection       # set the parameter file of an injection
    self._injection_parameters = None
    self._injection_credible_regions = {} # dictionary of minimal credible regions within which an injected parameter is found
    self._harmonics = harmonics       # list of frequency harmonics
    self._modeltype = modeltype       # the model type ('waveform' or 'source')
    self._biaxial = biaxial           # whether the source is a biaxial star (rather than triaxial)
    self._usegwphase = usegwphase     # whether to use GW phase rather than rotational phase

    # check if injection parameter file has been given
    if self._injection != None:
      # try and read it
      try:
        self._injection_parameters = psr_par(self._injection)
      except:
        print("Error... cannot read injection parameter file '%s'." % self._injection, file=sys.stderr)
        sys.exit(1)

      if self._usegwphase: # change initial phase if required
        if hasattr(self._injection_parameters, 'PHI0'):
          phi0val = 2.*self._injection_parameters['PHI0']
          setattr(self._injection_parameters, 'PHI0', phi0val)

    if 'Joint' in self._ifos: # put 'Joint' at the end
      j = self._ifos.pop(self._ifos.index('Joint'))
      self._ifos.append(j)

    # check files exist
    for ifo in postfiles:
      if not os.path.isfile(postfiles[ifo]):
        print("Error... posterior samples file '%s' for '%s' does not exist." % (postfiles[ifo], ifo), file=sys.stderr)
        sys.exit(1)
      else: # convert into posterior class
        pos, sigev, noiseev = pulsar_nest_to_posterior(postfiles[ifo])
        self._posteriors[ifo] = pos
        self._signal_evidence[ifo] = sigev
        self._noise_evidence[ifo] = noiseev
        self._Bsn[ifo] = sigev - noiseev

        theseparams = []
        for param in pos.names:
          thispar = True
          # remove parameters with 'logl', 'logw' or 'logprior' in the name
          for testpar in ['logl', 'logw', 'logprior']:
            if testpar in param.lower():
              thispar = False
              break
          if thispar:
            theseparams.append(param)

        if len(self._parameters) == 0:
          self._parameters = copy.copy(theseparams)
        else: # make sure different detectors do not have different sets of parameters
          for param in theseparams:
            if param not in self._parameters:
              print("Error... parameter '%s' is not defined in posteriors samples for '%s'." % (param, ifo), file=sys.stderr)
              sys.exit(1)

        if self._usegwphase: # try switching phi0 to 2*phi0 if working with l=m=2 gravitational wave initial phase (e.g. for hardware injections)
          if 'phi0' in pos.names:
            phi0new = bppu.PosteriorOneDPDF('phi0', 2.*pos['phi0'].samples)
            pos.pop('phi0')
            pos.append(phi0new)

        # if just working with results at 2f, but C22 and phi22 parameters have been used, convert back to h0 and phi0
        if len(self._harmonics) == 1:
          if self._harmonics[0] == 2.:
            if 'c22' in pos.names and self._modeltype == 'waveform':
              h0new = bppu.PosteriorOneDPDF('h0', 2.*pos['c22'].samples)
              pos.pop('c22')
              pos.append(h0new)
              ih0 = self._parameters.index('h0')
              self._parameters.pop(ih0)
              self._parameters.insert(ih0, 'C22')
              if self._injection != None:
                if not hasattr(self._injection_parameters, 'H0'):
                  setattr(self._injection_parameters, 'H0', 2.*self._injection_parameters['C22'])

            if 'phi22' in pos.names and modeltype == 'waveform':
              phi0new = bppu.PosteriorOneDPDF('phi0', 0.5*pos['phi22'].samples)
              pos.pop('phi0')
              pos.append(phi0new)
              iphi0 = self._parameters.index('phi0')
              self._parameters.pop(iphi0)
              self._parameters.insert(iphi0, 'PHI0')
              if self._injection != None:
                if not hasattr(self._injection_parameters, 'PHI0'):
                  setattr(self._injection_parameters, 'PHI0', 0.5*self._injection_parameters['PHI22'])

        # for a biaxial star (and if using the waveform model) set phi22 = 2*phi21
        if len(self._harmonics) == 2:
          if 1. in self._harmonics and 2. in self._harmonics and self._biaxial and self._modeltype == 'waveform':
            if 'phi21' in pos.names:
              phi22 = bppu.PosteriorOneDPDF('phi22', 2.*pos['phi21'].samples)
              pos.append(phi22)
              if self._injection != None:
                if hasattr(self._injection_parameters, 'PHI21'):
                  setattr(self._injection_parameters, 'PHI22', 2.*self._injection_parameters['PHI21'])

        self._posteriors[ifo] = pos

    self._get_snrs() # get the optimal SNR of the signal

    self._get_bayes_factors() # get the Bayes factors

  @property
  def parameters(self):
    # return a list of the source parameter names
    return self._parameters

  @property
  def bsn(self):
    # return a dictionary of the signal vs noise log Bayes factors
    return self._Bsn

  @property
  def snrs(self):
    # return dictionary of snrs
    return self._optimal_snrs

  @property
  def bci(self):
    # return coherent vs incoherent Bayes factor
    return self._Bci

  @property
  def bcin(self):
    # return coherent vs incoherent or noise Bayes factor
    return self._Bcin

  def h0_ul(self, ifo):
    # return the h0 upper limit for a given detector
    if ifo in self._h0ul:
      return self._h0ul[ifo]
    else:
      return None

  def ellipticity_ul(self, ifo):
    # return the ellipticity upper limit
    if ifo in self._ellipticity:
      return self._ellipticity[ifo]
    else:
      return None

  def q22_ul(self, ifo):
    # return the Q22 quadrupole moment upper limit
    if ifo in self._q22:
      return self._q22[ifo]
    else:
      return None

  def C21_ul(self, ifo):
    # return the C21 upper limit
    if ifo in self._C21:
      return self._C21[ifo]
    else:
      return None

  def C22_ul(self, ifo):
    # return the C22 upper limit
    if ifo in self._C22:
      return self._C22[ifo]
    else:
      return None

  def I21_ul(self, ifo):
    # return the I21 upper limit
    if ifo in self._I21:
      return self._I21[ifo]
    else:
      return None

  def I31_ul(self, ifo):
    # return the I31 upper limit
    if ifo in self._I31:
      return self._I31[ifo]
    else:
      return None

  def sdlim_ratio(self, ifo):
    # return the spin-down limit ratio
    if ifo in self._sdratio:
      return self._sdratio[ifo]
    else:
      return None

  def _get_snrs(self):
    # get the returned match filter SNRs
    for ifo in self._ifos:
      self._optimal_snrs[ifo] = self.get_snr(os.path.dirname(self._postfiles[ifo]))

  def get_snr(self, pdir):
    # get SNR from files in pdir
    snr = 0.
    # get files with SNR in the name in the posterior file directory
    snrfiles = [sf for sf in os.listdir(pdir) if 'SNR' in sf]
    if len(snrfiles) < 1:
      print("Error... no SNR files are given for '%s'" % ifo, file=sys.stderr)
      sys.exit(1)

    for snrfile in snrfiles: # average SNR values from mulitple files (e.g. if nested sampling has been run multiple times for the same signal)
      fp = open(os.path.join(pdir, snrfile), 'r')
      lines = [line.strip() for line in fp.readlines()]

      if '# Recovered SNR' not in lines:
        print("Error... no recovered SNRs are given in the SNR file '%s'." % snrfile, file=sys.stderr)
        sys.exit(1)

      # just return the final (either coherent or single detector-single frequency value of the SNR)
      linevals = lines[-1].split()
      snr += float(linevals[-1])
    return snr/len(snrfiles)

  def _get_bayes_factors(self):
    # get the Bayes factors for the signal
    incoherent = 0. # the incoherent evidence
    for ifo in self._ifos:
      self._Bsn[ifo], self._signal_evidence[ifo], self._noise_evidence[ifo], self._maxL[ifo] = self.get_bayes_factor(self._postfiles[ifo])

      if ifo != 'Joint':
        incoherent += self._signal_evidence[ifo]

    # get the coherent vs incoherent noise evidence
    if len(self._ifos) > 2 and 'Joint' in self._ifos:
      self._Bci = self._signal_evidence['Joint'] - incoherent
      self._Bcin = self._signal_evidence['Joint'] - np.logaddexp(incoherent, self._noise_evidence['Joint'])

  def get_bayes_factor(self, postfile):
    # return the Bayes factor extracted from a posterior file
    try:
      fe = os.path.splitext(postfile)[-1].lower() # file extension
      if fe == '.h5' or fe == '.hdf': # HDF5 file
        # open hdf5 file
        f = h5py.File(postfile, 'r')
        a = f['lalinference']['lalinference_nest']
        evdata = (a.attrs['log_bayes_factor'], a.attrs['log_evidence'], a.attrs['log_noise_evidence'], a.attrs['log_max_likelihood'])
        f.close()
      else: # try old/legacy file format
        B = np.loadtxt(postfile.replace('.gz', '')+'_B.txt')
        evdata = tuple(B.tolist())
    except:
      print("Error... could not extract evidences from '%s'." % postfile, file=sys.stderr)
      sys.exit(1)

    return evdata # line contains Bsn, signal evidence, noise evidence, max likelihood

  def snr(self, ifo):
    # return the SNR for a given detector
    if ifo in self._optimal_snrs:
      return self._optimal_snrs[ifo]
    else:
      return None

  def create_joint_posterior_plot(self, parameters, bins=20, ifo=None, truths=None, credintervals=[0.9], filename=None,
                                  figformats=['png'], ratio=3, figlimits=None, contourlimits=None, jointsamples=True, whichtruth=None, scatter_kwargs={}):
    # create a plot with the 1D and 2D joint posteriors (if ifo is None then use all ifos in class)
    if ifo != None and ifo in self._ifos:
      plotifos = [ifo]
    else:
      plotifos = self._ifos

    # if a joint detector posterior is present, make that first so other data is plotted on top
    if 'Joint' in plotifos:
      j = plotifos.pop(plotifos.index('Joint'))
      plotifos.insert(0, j)

    coldict = {'H1': 'red', 'H2': 'cyan', 'L1': 'green', 'V1': 'blue', 'G1': 'magenta', 'Joint': 'black'}

    if len(parameters) < 2 or len(parameters) > len(self._parameters):
      print("Error... can only plot posterior distributions for more than one parameters", file=sys.stderr)
      sys.exit(1)

    thistruth = []
    for param in parameters:
      if param not in self._parameters:
        print("Error... the requested parameter '%s' is not recognised" % param, file=sys.stderr)
        sys.exit(1)

      if self._injection_parameters != None and truths == None: # we have injection values
        thistruth.append(self._injection_parameters[param.upper()])

    if self._injection_parameters != None and truths == None: truths = thistruth

    if truths != None:
      if isinstance(truths, list): # just a list of the true parameter values
        newtruths = {}
        if len(truths) != len(parameters):
          print("Error... number of true values must be the same as number of parameters.", file=sys.stderr)
          sys.exit(1)
        else:
          # convert to dictionary
          for ifo in plotifos:
            newtruths[ifo] = truths
          truths = newtruths
      if isinstance(truths, dict): # dictionary keyed to IFO with lists of true parameters for each IFO (this is really designed for the SNR vs Bayes factor plots)
        for ifo in plotifos:
          if ifo not in truths:
            print("Error... problem in truths values. Detector '%s' not recognised." % ifo, file=sys.stderr)
            sys.exit(1)
          else:
            if len(truths[ifo]) != len(parameters):
              print("Error... number of true values must be the same as number of parameters.", file=sys.stderr)
              sys.exit(1)
    else:
      truths = {}
      for ifo in plotifos:
        truths[ifo] = None # set all to be None

    # use first ifo and get the required posterior samples
    x = self._posteriors[plotifos[0]][parameters[0]].samples
    labels = []
    if parameters[0].upper() in paramlatexdict:
      labels.append(paramlatexdict[parameters[0].upper()])
    else:
      labels.append(parameters[0])
    for param in parameters[1:]:
      x = np.hstack((x, self._posteriors[plotifos[0]][param].samples))
      if param.upper() in paramlatexdict:
        labels.append(paramlatexdict[param.upper()])
      else:
        labels.append(param)

    # set styles to different detectors
    if plotifos[0] == 'Joint':
      histops = {'histtype': 'stepfilled', 'color': 'darkslategrey', 'edgecolor': coldict['Joint'], 'linewidth': 1.5}
      contourops = {'colors': coldict[plotifos[0]]}
      showcontours = True
      if whichtruth == 'Joint' or whichtruth == 'all':
        truthops = {'color': 'black', 'markeredgewidth': 2}
      else:
        truthops = {}
      showpoints = jointsamples
    else:
      showcontours = True
      contourops = {'colors': 'dark'+coldict[plotifos[0]]}
      showpoints = True
      if len(plotifos) == 1: # if just one detector use a filled histogram
        histops = {'histtype': 'stepfilled', 'color': coldict[plotifos[0]], 'edgecolor': coldict[plotifos[0]], 'linewidth': 1.5}
        truthops = {'color': 'black', 'markeredgewidth': 2}
      else:
        histops = {'histtype': 'step', 'color': coldict[plotifos[0]], 'linewidth': 1}
        if whichtruth == plotifos[0]:
          truthops = {'color': 'black', 'markeredgewidth': 2}
        elif whichtruth == 'all':
          truthops = {'color': 'dark'+coldict[plotifos[0]], 'markeredgewidth': 2}
        else:
          truthops = {}

    sc = scotchcorner(x, bins=bins, ratio=ratio, labels=labels, truths=truths[plotifos[0]], datatitle=plotifos[0], showlims='both', hist_kwargs=histops,
                      showcontours=showcontours, contour_levels=credintervals, contour_kwargs=contourops, truths_kwargs=truthops, contour_limits=contourlimits,
                      limits=figlimits, show_level_labels=False, showpoints=showpoints, scatter_kwargs=scatter_kwargs)

    # now add the rest to the plots
    if len(plotifos) > 1:
      for k, ifo in enumerate(plotifos[1:]):
        histops = {'histtype': 'step', 'color': coldict[ifo], 'linewidth': 1}
        x = self._posteriors[ifo][parameters[0]].samples
        for param in parameters[1:]:
          x = np.hstack((x, self._posteriors[ifo][param].samples))
        showcontours = True
        contourops = {'colors': 'dark'+coldict[plotifos[k+1]]}
        if whichtruth == plotifos[k+1]:
          truthops = {'color': 'black', 'markeredgewidth': 2}
        elif whichtruth == 'all':
          truthops = {'color': 'dark'+coldict[plotifos[k+1]], 'markeredgewidth': 2}
        else:
          truthops = {}
        sc.add_data(x, hist_kwargs=histops, datatitle=ifo, truths=truths[ifo], showcontours=showcontours, contour_kwargs=contourops, contour_levels=credintervals, show_level_labels=False, truths_kwargs=truthops, scatter_kwargs=scatter_kwargs, contour_limits=contourlimits, limits=figlimits)

    # output the plots
    if 'png' not in figformats and 'svg' not in figformats:
      print("Error... must include 'png' and/or 'svg' in the output figure formats", file=sys.stderr)
      sys.exit(1)

    if filename == None: # use default file names
      if len(parameters) == len(self._parameters):
        outfilepre = os.path.join(self._outputdir, 'all_posteriors')
      else:
        outfilepre = os.path.join(self._outputdir, "vs".join(parameters))
    else:
      outfilepre = os.path.join(self._outputdir, filename)

    outfiles = []
    for ftype in figformats:
      outfile = outfilepre+'.'+ftype
      try:
        sc.fig.subplots_adjust(left=0.18, bottom=0.15) # adjust size
        sc.savefig(outfile)
      except:
        print("Error... could not output posterior plot file '%s'." % outfile, file=sys.stderr)
        sys.exit(1)
      outfiles.append(outfile)

    return outfiles # list of output figure filenames

  # create standard joint plots
  def create_joint_plots_table(self, allparams=False, title='Joint distributions'):
    # create a table containing a set of standard 2D or 3D joint posterior plots (depending on the analysis setup)
    header = htmltag('h2', title, newline=True)

    table = htmltable()
    table.addrow()
    paramlist = []

    # set limits for credible interval contours plots using maximum allowed ranges (see e.g. Tables 1 and 2 of http://arxiv.org/abs/1508.00416v2)
    limits = {}
    for param in self._parameters:
      if param == 'h0':
        # only set lower limit on amplitudes if posteriors are close (less than 3 sigma) to zero
        limits[param] = ()
        for ifo in self._ifos:
          if self._posteriors[ifo][param].mean - 3.*self._posteriors[ifo][param].stdev < 0.:
            limits[param] = (0., None) # must be greater than 0.
            break
      elif param == 'cosiota':
        limits[param] = (-1., 1.) # cos(iota)
      elif param == 'phi0':
        if self._usegwphase:
          limits[param] = (0., 2.*np.pi)
        else:
          limits[param] = (0., np.pi)
      elif param == 'psi':
        if self._biaxial and self._modeltype == 'source':
          limits[param] = (0., np.pi)
        else:
          limits[param] = (0., np.pi/2.)
      elif param == 'c21' or param == 'c22':
        if self._biaxial and self._modeltype == 'waveform':
          limits[param] = ()
        else:
          limits[param] = ()
          for ifo in self._ifos:
           if self._posteriors[ifo][param].mean - 3.*self._posteriors[ifo][param].stdev < 0.:
             limits[param] = (0., None) # must be greater than 0.
             break
      elif param == 'lambda':
        limits[param] = (0., np.pi)
      elif param == 'costheta':
        limits[param] = (0., 1.)
      else: # default to empty tuple
        limits[param] = () # empty tuple

    if allparams: # plot all parameters
      paramlist = [self._parameters]
    else:
      if len(self._harmonics) == 1:
        if harmonics[0] == 2.: # for emission at twice the rotation frequency
          # plot h0 versus cos(iota) and psi verus phi0
          paramlist = [['h0', 'cosiota'], ['phi0', 'psi']]
        else:
          if self._modeltype == 'source':
            print("Error... do not know 'source' parameterisation for emission at purely the rotation frequency", file=sys.stderr)
            sys.exit(1)
          else:
            # plot C21 versus cosiota and psi versus phi21
            paramlist = [['c21', 'cosiota'], ['phi21', 'psi']]
      else:
        if self._modeltype == 'source':
          if self._biaxial:
            # plot I31 vs cosiota vs phi0 and psi vs lambda vs theta (I21 is zero)
            paramlist = [['i31', 'cosiota', 'phi0'], ['psi', 'lambda', 'costheta']]
          else:
            # plot I21 vs I31 vs cosiota and phi0 vs psi vs lambda vs theta
            paramlist = [['i21', 'i31', 'cosiota'], ['phi0', 'psi', 'lambda', 'costheta']]
        else:
          # plot C21 vs C22 vs cosiota and psi vs phi21 vs phi22
          paramlist = [['c21', 'c22', 'cosiota'], ['psi', 'phi21', 'phi22']]

    notpresent = False
    for parampairs in paramlist:
      contourlimits = []
      figlimits = []
      for p in parampairs:
        if p not in self._parameters:
          notpresent = True
        else:
          contourlimits.append(limits[p])
          figlimits.append(limits[p])
      if notpresent: break
      # contourlimits = None # temporary for testing
      pf = self.create_joint_posterior_plot(parampairs, figformats=['png'], ratio=2, figlimits=figlimits, contourlimits=contourlimits, jointsamples=False)
      if allparams:
        tagclass = 'jointplot'
      else:
        tagclass = 'posplot'
      table.adddata(atag(os.path.basename(pf[0]), '<img class="{}" src="{}"/>'.format(tagclass, os.path.basename(pf[0]))).text)

    return header.text + table.tabletext

  def create_sample_plot_table(self, figformats=['png']):
    # create the table of posterior sample plots
    if 'png' not in figformats:
      print("Error... must include 'png' in the output figure formats", file=sys.stderr)
      sys.exit(1)

    header = htmltag('h2', 'Posterior samples', newline=True) # the output html table

    table = htmltable()

    # get the plots
    for param in self._parameters:
      chainfig = plot_posterior_chain(self._posteriors.values(), param, self._ifos, grr=None)
      chainfile = os.path.join(self._outputdir, param+'_postchain')

      for ftype in figformats:
        thischainfile = chainfile+'.'+ftype
        try:
          chainfig.savefig(thischainfile)
        except:
          print("Error... could not output posterior samples plot for '%s'." % param, file=sys.stderr)
          sys.exit(1)

        if ftype == 'png':
          table.addrow()
          table.adddata(atag(os.path.basename(thischainfile), '<img class="chainplot" src="{}"/>'.format(os.path.basename(thischainfile))).text)

    self._samples_table = header.text + table.tabletext
    return self._samples_table # output html for table

  def create_stats_table(self, credints=[95]):
    # create table of posterior sample statistics, including credible intervals given by the list
    self._credint = credints
    for ci in self._credint:
      if ci >= 100. or ci <= 0.:
        print("Error... credible interval '%s' is outside the range between 0%% and 100%%." % str(ci), file=sys.stderr)
        sys.exit(1)

    header = htmltag('h2', 'Posterior statistics', newline=True)

    table = htmltable()

    cols = ['Detector', 'max. poterior', 'mean', 'median', '&sigma;']
    for ci in self._credint:
      cols.append('%d%% credible interval' % int(ci))
    if self._injection_parameters != None:
      cols.insert(0, 'Inj. value') # prepend the injection value
    cols.insert(0, '&nbsp;') # add empty cell as header line for the parameter name

    # create header row
    table.addrow()
    for col in cols:
      table.adddata(col, header=True)

    # loop through parameters header row to table
    for i, param in enumerate(self._parameters):
      pu = param.upper()
      dispkwargs = {} # any additional arguments required by value display function

      # parameter html to display
      if pu in paramhtmldict:
        pdisp = paramhtmldict[pu]

        # some different display styles (compared to the defaults) for certain parameters
        if pu == 'RA' or pu == 'DEC' or pu == 'RAJ' or pu == 'DECJ':
          dispkwargs = {'stype': 'rads'} # display as rad rather than hh/dd:mm:ss string
        if pu == 'F0': dispkwargs = {'dp': 2} # just display with 2 decimal places
      else:
        pdisp = param

      # parameter value function to display
      if pu in paramhtmldispfunc.__dict__:
        dispfunc = paramhtmldispfunc.__dict__[pu]
      else:
        dispfunc = paramhtmldispfunc.__dict__['DEFAULTSTR']

      for j, ifo in enumerate(self._ifos):
        table.addrow()
        if j == 0: # add parameter name and injection value first)
          table.adddata(pdisp, rowspan=len(self._ifos))
          if self._injection_parameters != None:
            if hasattr(self._injection_parameters, param.upper()):
              table.adddata(dispfunc(str(self._injection_parameters[param.upper()]), **dispkwargs), rowspan=len(self._ifos))
            else:
              table.adddata("*", rowspan=len(self._ifos)) # no value for this parameter

        if i == 0:
          self._injection_credible_regions[ifo] = {} # dictionary for injection minimal credible regions to be input for each parameter
        maxL, maxLparams = self._posteriors[ifo].maxL
        table.adddata(ifo, dataclass=ifo)
        table.adddata(dispfunc(str(maxLparams[param]), **dispkwargs))    # maximum likelihood
        table.adddata(dispfunc(str(self._posteriors[ifo].means[param]), **dispkwargs))   # mean value
        table.adddata(dispfunc(str(self._posteriors[ifo].medians[param]), **dispkwargs)) # median value
        table.adddata(dispfunc(str(self._posteriors[ifo].stdevs[param]), **dispkwargs))  # standard deviations

        for k, ci in enumerate(credints):
          paramval = None
          if k == 0 and self._injection_parameters != None:
            if hasattr(self._injection_parameters, param.upper()):
              paramval = self._injection_parameters[param.upper()]
          low, high, cr = self.credible_interval(ifo, param, ci, paramval)

          if k == 0 and self._injection_parameters != None:
            if hasattr(self._injection_parameters, param.upper()):
              self._injection_credible_regions[ifo][param] = cr

          table.adddata('({}, {})'.format(dispfunc(str(low), **dispkwargs), dispfunc(str(high), **dispkwargs)))

    self._stats_section = header.text + table.tabletext
    return self._stats_section

  def create_limits_table(self, freq, sdlim=None, dist=None, ul=95):
    # create a table with upper limits on amplitudes and Bayes factors
    self._limits_table = ""

    table = htmltable()

    self._h0ul = {}
    self._C21 = {}
    self._C22 = {}
    self._I21 = {}
    self._I31 = {}
    self._sdratio = {}
    self._ellipticity = {}
    self._q22 = {}

    # add header row
    table.addrow()

    table.adddata('&nbsp;', header=True) # nothing in the first column

    if len(self._harmonics) == 1 and self._harmonics[0] == 2:
      table.adddata(paramhtmldict['H0UL'].format(ul), header=True)
      table.adddata(paramhtmldict['ELL'], header=True)
      table.adddata(paramhtmldict['Q22'], header=True)
      table.adddata(paramhtmldict['SDRAT'], header=True)
    elif len(self._harmonics) == 1 and self._harmonics[0] == 1 and self._modeltype == 'waveform':
      table.adddata(paramhtmldict['C21UL'].format(ul), header=True)
    elif len(self._harmonics) == 2 and self._modeltype == 'waveform':
      table.adddata(paramhtmldict['C21UL'].format(ul), header=True)
      table.adddata(paramhtmldict['C22UL'].format(ul), header=True)
      table.adddata(paramhtmldict['H0UL'].format(ul) + ' (from C<sub>22</sub>)', header=True)
      table.adddata(paramhtmldict['ELL'] + ' (from C<sub>22</sub>)', header=True)
      table.adddata(paramhtmldict['Q22'] + ' (from C<sub>22</sub>)', header=True)
      table.adddata('ratio', header=True)
    elif len(self._harmonics) == 2 and self._modeltype == 'source':
      if not self._biaixal:
        table.adddata(paramhtmldict['I21UL'].format(ul), header=True) # not needed for a biaxial source
      table.adddata(paramhtmldict['I31UL'].format(ul), header=True)
    else:
      print("Error... do not know how to output results table for this situation.", file=sys.stderr)
      sys.exit(1)

    # add SNR
    table.adddata(paramhtmldict['SNR'], header=True)

    # add max log likelihood
    table.adddata(paramhtmldict['MAXL'], header=True)

    # add Bayes factor (signal vs noise)
    table.adddata(paramhtmldict['BSN'], header=True)

    # check whether to add coherent vs incoherent Bayes factor (and coherent vs (incoherent or noise) Bayes factor)
    if len(self._ifos) > 2 and 'Joint' in self._ifos:
      table.adddata(paramhtmldict['BCI'], header=True)
      table.adddata(paramhtmldict['BCIN'], header=True)

    for ifo in self._ifos:
      table.addrow()
      table.adddata('{}'.format(ifo), dataclass=ifo)

      if self._modeltype == 'waveform' and len(self._harmonics) == 2:
        for p in ['c21', 'c22']:
          Cul = upper_limit_greedy(self._posteriors[ifo][p].samples, upperlimit=(ul/100.0))
          table.adddata('{}'.format(exp_str(Cul)), dataclass=ifo)
          if p == 'c21':
            self._C21[ifo] = Cul
          if p == 'c22':
            self._C22[ifo] = Cul

        self._h0ul[ifo] = self._C22[ifo]*2.
        table.adddata('{}'.format(exp_str(self._h0ul[ifo])), dataclass=ifo) # convert to h0

      if (len(self._harmonics) == 1 and self._harmonics[0] == 2) or (len(self._harmonics) == 2 and self._modeltype == 'waveform'):
        # get h0 upper limit
        if len(self._harmonics) == 1 and self._harmonics[0] == 2:
          h0ul = upper_limit_greedy(self._posteriors[ifo]['h0'].samples, upperlimit=(ul/100.0))
          self._h0ul[ifo] = h0ul
          table.adddata('{}'.format(exp_str(h0ul)), dataclass=ifo)

        # ellipticity
        if dist == None:
          self._ellipticity[ifo] = None
          table.adddata('*', dataclass=ifo)
        else:
          self._ellipticity[ifo] = h0_to_ellipticity(self._h0ul[ifo], freq, dist)
          table.adddata('{}'.format(exp_str(self._ellipticity[ifo])), dataclass=ifo)

        # quadrupole
        if dist == None:
          self._q22[ifo] = None
          table.adddata('*', dataclass=ifo)
        else:
          self._q22[ifo] = h0_to_quadrupole(self._h0ul[ifo], freq, dist)
          table.adddata('{}'.format(exp_str(self._q22[ifo])), dataclass=ifo)

         # get spin down limit ratio
        if sdlim == None:
          self._sdratio[ifo] = None
          table.adddata('*', dataclass=ifo)
        else:
          self._sdratio[ifo] = self._h0ul[ifo]/sdlim
          table.adddata('%.02f' % self._sdratio[ifo], dataclass=ifo)

      if len(self._harmonics) == 2 and self._modeltype == 'source':
        if not self._biaxial:
          self._I21[ifo] = upper_limit_greedy(self._posteriors[ifo]['I21'].samples, upperlimit=(ul/100.0))
          table.adddata('{}'.format(exp_str(self._I21[ifo])), dataclass=ifo)

        self._I31[ifo] = upper_limit_greedy(self._posteriors[ifo]['I31'].samples, upperlimit=(ul/100.0))
        table.adddata('{}'.format(exp_str(self._I31[ifo])), dataclass=ifo)

      # add SNR
      if ifo in self._optimal_snrs:
        table.adddata('%.1f' % self._optimal_snrs[ifo], dataclass=ifo)
      else:
        table.adddata('*', dataclass=ifo)

      # add max log likelihood values
      table.adddata('%.f' % (self._maxL[ifo]/np.log(10.)), dataclass=ifo) # convert to log base 10
      table.adddata('%.1f' % (self._Bsn[ifo]/np.log(10.)), dataclass=ifo) # convert to log base 10

      # add coherent vs incoherent and coherent vs (incoherent or noise) for Joint ifo)
      if len(self._ifos) > 2 and 'Joint' in self._ifos:
        if ifo == 'Joint':
          table.adddata('%.1f' % (self._Bci/np.log(10.)), dataclass=ifo)  # convert to log base 10
          table.adddata('%.1f' % (self._Bcin/np.log(10.)), dataclass=ifo) # convert to log base 10
        else:
          table.adddata('*', dataclass=ifo)
          table.adddata('*', dataclass=ifo)

    self._limits_table += table.tabletext
    return self._limits_table

  def credible_interval(self, ifo, param, ci=95, paramval=None):
    # get the given percentage (minimum) credible interval from samples for detector given by ifo and param
    # (for injections, where we have parameter values (paramval) get the corresponding smallest credible
    # region that contains the parameter value)
    samples = self._posteriors[ifo][param].samples.squeeze()
    samples.sort()
    lowbound, highbound, _ = self._ci_loop(samples, ci)

    cival = None
    if paramval != None:
      cifound = False
      # loop over different credible intervals until finding the one that contains the injection
      for cit in range(1, 101):
        l, h, cr = self._ci_loop(samples, cit)
        if paramval >= l and paramval <= h:
          cifound = True
          cival = cit
          break
      if not cifound:
        cival = 100

    return lowbound, highbound, cival

  def _ci_loop(self, sortedsamples, ci):
    lowbound = sortedsamples[0]
    highbound = sortedsamples[-1]
    cregion = highbound - lowbound
    lsam = len(sortedsamples)
    cllen = int(lsam*float(ci)/100.)
    for j in range(lsam-cllen):
      if sortedsamples[j+cllen] - sortedsamples[j] < cregion:
        lowbound = sortedsamples[j]
        highbound = sortedsamples[j+cllen]
        cregion = highbound - lowbound

    return lowbound, highbound, cregion


class create_background(posteriors):
  """
  Get information (evidence ratios and SNRs) from any the background analyses
  """
  def __init__(self, backgrounddirs, snrs, Bsn, outputdir, Bci=None, Bcin=None):
    # initialise with a dictionary (keyed to detectors) of directories containing the background analyses,
    # a dictionary of signal vs noise Bayes factors, a coherent vs incoherent Bayes factor and a coherent
    # vs incoherent or noise Bayes factor (all created from the posterior class)
    self._ifos = list(backgrounddirs.keys()) # detectors
    self._backgrounddirs = backgrounddirs
    self._dir_lists = {}    # dictionary with a list of background results directories for each detector
    self._Bci_fore = Bci    # the "foreground" Bayes factor for coherent vs incoherent
    self._Bcin_fore = Bcin  # the "foreground" Bayes factor for coherent vs incoherent or noise

    self._Bsn = {}          # dictionary of "background" signal vs noise Bayes factor lists for each detector
    self._Bci = []          # list of "background" coherent vs incoherent Bayes factors
    self._Bcin = []         # list of "background" coherent vs incoherent or noise Bayes factors
    self._optimal_snrs = {} # dictionary of "background" SNR lists for each detector
    self._signal_evidence = {}
    self._noise_evidence = {}
    self._maxL = {}

    # figure and contour limits
    self._figlimits = [(0., None), ()] # assumes SNR is first value
    self._contourlimits = [(0., None), ()]

    # set some default arguments for create_joint_posterior_plot
    self._plot_kwargs = {'ratio': 2, 'whichtruth': 'all', 'scatter_kwargs': {'alpha': 0.5}, 'figlimits': self._figlimits, 'contourlimits': self._contourlimits}

    self._Bsn_prob = {}    # probability of the "foreground value" given a Gaussian KDE applied to the background
    self._Bci_prob = None
    self._Bcin_prob = None

    self._injection_parameters = None # set to None

    self._outputdir = outputdir # set output directory

    dirlen = []
    for i, ifo in enumerate(self._ifos):
      if ifo not in Bsn or ifo not in snrs:
        print("Error... Bsn/SNR dictionary is not consistent with the background directories dictionary.", file=sys.stderr)
        sys.exit(1)

      self._Bsn[ifo] = [] # create empty list
      self._optimal_snrs[ifo] = [] # create empty list
      self._signal_evidence[ifo] = []
      self._noise_evidence[ifo] = []
      self._maxL[ifo] = []

      # get directory lists
      self._dir_lists[ifo] = [os.path.join(self._backgrounddirs[ifo], d) for d in os.listdir(self._backgrounddirs[ifo]) if os.path.isdir(os.path.join(self._backgrounddirs[ifo], d))]
      dirlen.append(len(self._dir_lists[ifo]))
      if dirlen[i] == 0:
        print("Error... no background results directories were present for '%s'." % ifo, file=sys.stderr)
        sys.exit(1)

      # check there are the same number of background runs for each ifo
      if i > 0:
        if dirlen[i] != dirlen[0]:
          print("Error... number of background results directories not then same for different detectors." % ifo, file=sys.stderr)
          sys.exit(1)

    self._Bsn_fore = Bsn    # the "foreground" Bayes factor for signal vs noise
    self._snrs_fore = snrs  # the "foreground" optimal snrs

    # get SNRs
    self._get_snrs()

    # get Bayes factors
    self._get_bayes_factors()

    # get probabilities the background distribution being greater than the foreground value
    self._get_bsn_prob()
    if len(self._ifos) > 2 and 'Joint' in self._ifos:
      self._get_bci_prob()

  @property
  def bci_prob(self):
    return self._Bci_prob

  @property
  def bcin_prob(self):
    return self._Bcin_prob

  def bsn_prob(self, ifo):
    return self._Bsn_prob[ifo]

  def _get_snrs(self):
    # get the returned match filter SNRs
    for ifo in self._ifos: # loop over detectors
      # loop over background directories
      for pdir in self._dir_lists[ifo]:
        self._optimal_snrs[ifo].append(self.get_snr(pdir))

  def _get_bayes_factors(self):
    # get the Bayes factors for the signal
    incoherent = np.zeros(len(self._dir_lists[self._ifos[0]])) # the incoherent evidence

    for ifo in self._ifos:
      for i, pdir in enumerate(self._dir_lists[ifo]):
        pfiles = os.listdir(pdir)
        Bsn = None
        for pfile in pfiles: # get HDF5 file with extenstion '.hdf' or '.h5'
          if fnmatch.fnmatch(pfile, 'posterior_samples*.hdf') or fnmatch.fnmatch(pfile, 'posterior_samples*.h5'):
            Bsn, sigev, noiseev, maxL = self.get_bayes_factor(os.path.join(pdir, pfile))

            self._Bsn[ifo].append(Bsn)
            self._signal_evidence[ifo].append(sigev)
            self._noise_evidence[ifo].append(noiseev)
            self._maxL[ifo].append(maxL)
            break
        if Bsn == None:
          print("Error... no HDF5 file with 'posterior_samples' in the name was found in '%s'." % pdir, file=sys.stderr)
          sys.exit(1)

        if ifo != 'Joint':
          incoherent[i] += self._signal_evidence[ifo][i]

    # get the coherent vs incoherent noise evidence
    if len(self._ifos) > 2 and 'Joint' in self._ifos:
      for i in range(len(self._dir_lists[self._ifos[0]])):
        self._Bci.append(self._signal_evidence['Joint'][i] - incoherent[i])
        self._Bcin.append(self._signal_evidence['Joint'][i] - np.logaddexp(incoherent[i], self._noise_evidence['Joint'][i]))

  def _get_bsn_prob(self):
    # get the probability of the background distribution of signal vs noise Bayes factors
    # being greater than the foreground value for each detector by using a Gaussian KDE on the
    # 1D samples from the background
    for ifo in self._ifos:
      kernel = stats.gaussian_kde(self._Bsn[ifo])
      self._Bsn_prob[ifo] = kernel.integrate_box_1d(self._Bsn_fore[ifo], np.inf) # integrate between foreground value and infinity

  def _get_bci_prob(self):
    kernel = stats.gaussian_kde(self._Bci)
    self._Bci_prob = kernel.integrate_box_1d(self._Bci_fore, np.inf)
    kernel = stats.gaussian_kde(self._Bcin)
    self._Bcin_prob = kernel.integrate_box_1d(self._Bcin_fore, np.inf)

  def bsn_plot(self, credint=[0.5, 0.95]):
    # create 2D scatter plots of SNR vs Bsn (with KDE estimated probability contours - if present just for Joint)
    # and 1D histograms for each - use the create_joint_posterior_plot's function

    self._posteriors = {}
    truths = {}

    # create fake posterior objects (requires logL to be there) for each detector
    for ifo in self._ifos:
      self._posteriors[ifo] = bppu.Posterior((['tmp', 'logL'], np.zeros((100,2))))

      # add the SNR and Bsn values
      snrs = bppu.PosteriorOneDPDF('snr', np.array([self._optimal_snrs[ifo]]).T)
      self._posteriors[ifo].append(snrs)

      bsn = bppu.PosteriorOneDPDF('bsn', np.array([self._Bsn[ifo]]).T/np.log(10.)) # convert into base 10 log
      self._posteriors[ifo].append(bsn)

      # remove temporary varaibles
      self._posteriors[ifo].pop('tmp')
      self._posteriors[ifo].pop('logl')

      self._parameters = ['snr', 'bsn'] #  set parameters
      truths[ifo] = [self._snrs_fore[ifo], self._Bsn_fore[ifo]/np.log(10.)]

    return self.create_joint_posterior_plot(['snr', 'bsn'], bins=int(np.log2(len(snrs.samples))), truths=truths, credintervals=credint, filename='bsn', **self._plot_kwargs)

  def bci_plot(self, credint=[0.5, 0.95], which='bci'):
    self._posteriors = {}
    if which == 'bcin':
      truths = [self._snrs_fore['Joint'], self._Bcin_fore/np.log(10.)]
      self._parameters = ['snr', 'bcin']
    else: # default to bci
      truths = [self._snrs_fore['Joint'], self._Bci_fore/np.log(10.)]
      self._parameters = ['snr', 'bci']

      if which != 'bci': # force which to be 'bci' (as there are only the two options)
        which = 'bci'

    self._posteriors['Joint'] = bppu.Posterior((['tmp', 'logL'], np.zeros((100,2))))

    # add the SNR and Bci values
    snrs = bppu.PosteriorOneDPDF('snr', np.array([self._optimal_snrs['Joint']]).T)
    self._posteriors['Joint'].append(snrs)

    if which == 'bcin':
      bci = bppu.PosteriorOneDPDF('bcin', np.array([self._Bcin]).T/np.log(10.))
    else:
      bci = bppu.PosteriorOneDPDF('bci', np.array([self._Bci]).T/np.log(10.))
    self._posteriors['Joint'].append(bci)

    curifos = self._ifos # save current detectors
    self._ifos = ['Joint'] # set for just the joint detector
    bciplot = self.create_joint_posterior_plot(self._parameters, bins=int(np.log2(len(snrs.samples))), truths=truths, credintervals=credint, filename=which, **self._plot_kwargs)
    self._ifos = curifos # reset detectors (in case they are needed later)

    return bciplot

  def create_background_table(self):
    # create table with background plots in them
    background_plots_table = htmltag('h2', 'Evidence distributions', newline=True).text

    table = htmltable()
    table.addrow()
    cols = 1 # number of columns

    # add Bsn plot
    bsnplot = self.bsn_plot()
    table.adddata(atag(os.path.basename(bsnplot[0]), '<img class="backgroundplot" src="{}"/>'.format(os.path.basename(bsnplot[0]))).text)

    if self._Bci_fore != None:
      bciplot = self.bci_plot(which='bci')
      table.adddata(atag(os.path.basename(bciplot[0]), '<img class="backgroundplot" src="{}"/>'.format(os.path.basename(bciplot[0]))).text)
      cols += 1

    if self._Bcin_fore != None:
      bcinplot = self.bci_plot(which='bcin')
      table.adddata(atag(os.path.basename(bcinplot[0]), '<img class="backgroundplot" src="{}"/>'.format(os.path.basename(bcinplot[0]))).text)
      cols += 1

    # add probabilities for the background distribution being greater than the foreground value
    innertable = htmltable(tableclass="evidencetable")
    innertable.addrow()
    innertable.adddata('&nbsp;')
    innertable.adddata('Probability of background distribution being greater than foreground', colspan=(len(self._ifos)+cols-1), header=True)
    innertable.addrow()
    innertable.adddata('Detector', header=True)
    for ifo in self._ifos:
      innertable.adddata(ifo, dataclass=ifo)
    if self._Bci_fore != None:
      innertable.adddata('B<sub>CvI</sub>')
    if self._Bcin_fore != None:
      innertable.adddata('B<sub>CvIN</sub>')

    innertable.addrow()
    innertable.adddata('&nbsp;')
    for ifo in self._ifos:
      innertable.adddata(exp_str(self.bsn_prob(ifo)), dataclass=ifo)

    if self._Bci_fore != None:
      innertable.adddata(exp_str(float(self.bci_prob), 1))

    if self._Bcin_fore != None:
      innertable.adddata(exp_str(float(self.bcin_prob), 1))

    table.addrow()
    table.adddata(innertable.tabletext, colspan=cols)

    background_plots_table += table.tabletext

    return background_plots_table


if __name__=='__main__':
  description = """This script will create a results page for a single pulsar from the known pulsar analysis pipeline. A configuration (.ini) file is required."""
  epilog = """An example configuration file could contain the following:

# a section for general analysis information
[general]
parfile = 'path_to_par_file'  # the path to the pulsar parameter file
detectors = ['H1', 'L1']      # a list of the detectors to use
with_joint = True             # a boolean stating whether to also add the joint multi-detector analysis
joint_only = False            # a boolean stating whether to only output the joint multi-detector analysis
with_background = False       # a boolean stating whether to include background analyses
injection = False             # a boolean stating whether this pulsar is a software/hardware injection
upper_limit = 95              # the percentage credible upper limit value
credible_interval = [95]      # a list of the percentage credible intervals for output statistics
use_gw_phase = False          # a boolean stating whether to assume the initial phase parameter is the rotational, or gravitational wave (l=m=2), phase (e.g. if looking a hardware injections)
harmonics = [2]               # a list of the frequency harmonics used in this analysis
model_type = waveform         # either 'waveform' (default) or 'source' specify the parameterisation
biaxial = False               # set whether the signal searched for was from a biaxial source

# a section for parameter estimation files
[parameter_estimation]
posteriors = {'H1': 'path_to_H1_posteriors', 'L1': 'path_to_L1_posteriors', 'Joint': 'path_to_Joint_posteriors'}            # a dictionary (keyed on detectors) pointing to the locations of posterior sample files
background = {'H1': 'path_to_H1_background_dir', 'L1': 'path_to_L1_background_dir', 'Joint': 'path_to_Joint_backgroud_dir'} # a dictionary (keyed on detectors) pointing to directories containing the background posterior files

# a section for pre-processed data information
[data]
files = {'H1': 'path_to_H1_data', 'L1': 'path_to_L1_data'} # a dictionary (keyed on detectors) pointing to the locations (or lists of locations for multiple harmonics) of pre-processed data files

# a section for the output location for this pulsar
[output]
path = 'path_to_output_base_directory' # the path to the base directory in which the results page will be created
indexpage = 'path_to_index_page'       # an optional path (relative to the base directory) to the index page containing results from multiple pulsars

# a section for plotting options
[plotting]
all_posteriors = False # a boolean stating whether to show joint posterior plots of all parameters (default: False)
eps_output = False     # a boolean stating whether to also output eps versions of figures (png figures will automatically be produced)
pdf_output = False     # a boolean stating whether to also output pdf versions of figures (png figures will automatically be produced)

"""

  parser = argparse.ArgumentParser( description=description, epilog=epilog, formatter_class=argparse.RawDescriptionHelpFormatter )
  parser.add_argument("inifile", help="The configuration (.ini) file")

  # parse input options
  opts = parser.parse_args()

  inifile = opts.inifile

  # open and parse config file
  cp = ConfigParser()
  try:
    cp.read(inifile)
  except:
    print("Error... cannot parse configuration file '%s'" % inifile, file=sys.stderr)
    sys.exit(1)

  # get the output directory
  try:
    outdir = cp.get('output', 'path')
  except:
    print("Error... no output directory 'path' specified in [output] section.", file=sys.stderr)
    sys.exit(1)

  # create directory if required
  if not os.access(outdir, os.W_OK) and not os.path.isdir(outdir): # check if directory already exists
    try:
      os.makedirs(outdir)
    except:
      print("Error... cannot make output directory '%s'." % outdir, file=sys.stderr)
      sys.exit(1)

  # get the index page location
  try:
    indexpage = cp.get('output', 'indexpage')
  except:
    indexpage = None

  # get the pulsar parameter file
  parfile = None
  try:
    parfile = cp.get('general', 'parfile')
  except:
    print("Error... Pulsar parameter 'parfile' must specified in the [general] section.", file=sys.stderr)
    sys.exit(1)

  # read in data from par file
  try:
    par = psr_par(parfile)
  except:
    print("Error... Pulsar parameter (.par) file '%s' could not be opened!" % parfile, file=sys.stderr)
    sys.exit(1)

  # get pulsar PSRJ name
  pname = par['PSRJ']
  if not pname:
    print("Error... no PSRJ name set in pulsar parameter file '%s'." % parfile, file=sys.stderr)

  # check/get required parameters
  f0 = par['F0'] # pulsar rotation frequency (Hz)
  if not f0:
    print("Error... no 'F0' value in the parameter file '%s'." % parfile, file=sys.stderr)
    sys.exit(1)

  f1 = par['F1'] # pulsar first time derivative of rotational frequency (Hz/s)
  if not f1:
    f1 = 0. # set to zero if not given

  # get the upper limit credible interval
  try:
    upperlim = ast.literval_eval(cp.get('general', 'upper_limit'))
  except: # default to 95%
    upperlim = 95

  # get credible intervals for output statistics
  try:
    credints = ast.literval_eval(cp.get('general', 'credible_interval'))
  except: # default to 95%
    credints = [95]

  # check whether looking at an injection or not
  try:
    injection = cp.getboolean('general', 'injection')
  except:
    injection = False # if nothing is given then assume that this is not an injection

  if injection:
    injectionfile = parfile
  else:
    injectionfile = None

  # check whether to use the rotational, or gravitational wave phase (e.g. for hardware injections), in output plots
  try:
    usegwphase = cp.getboolean('general', 'use_gw_phase')
  except:
    usegwphase = False

  # attempt to get pulsar distance, proper motion corrected age and any association (e.g. GC from the ATNF catalogue)
  dist = age = assoc = sdlim = f1sd = None
  atnfurl = None
  if not injection:
    pinfo = get_atnf_info(pname)
    if pinfo != None:
      dist, age, assoc = pinfo # unpack values
      atnfversion = '1.54'
      atnfurl = 'http://www.atnf.csiro.au/people/pulsar/psrcat/proc_form.php?version=' + atnfversion
      atnfurl += '&startUserDefined=true&pulsar_names=' + re.sub('\+', '%2B', pname)
      atnfurl += '&ephemeris=long&submit_ephemeris=Get+Ephemeris&state=query'

    # if distance is in the par file use that instead
    if par['DIST']:
      dist = par['DIST']

    # set the corrected spin-down value (based on intrinsic age or using consverative value from GC pulsars)
    f1sd = set_spin_down(age, assoc, f0, f1)

    # get spin-down limit
    if f1sd != None and dist != None:
      sdlim = spin_down_limit(f0, f1sd, dist)

  # check whether to only include results from a joint (multi-detector) analysis
  try:
    jointonly = cp.getboolean('general', 'joint_only')
  except:
    jointonly = False

  # check which frequency harmonics were used in this analysis
  try:
    harmonics = ast.literal_eval(cp.get('general', 'harmonics'))
  except:
    harmonics = [2] # default to twice the rotation frequency

  # check whether the waveform or source model is being used
  try:
    modeltype = cp.get('general', 'model_type')
  except:
    modeltype = 'waveform'

  if modeltype not in ['waveform', 'source']:
    print("Error... unknown 'model type' '%s' specified." % modeltype, file=sys.stderr)
    sys.exit(1)

  try:
    biaxial = cp.getboolean('general', 'biaxial')
  except:
    biaxial = False

  # check whether a background analysis has also been performed
  try:
    with_background = cp.getboolean('general', 'with_background')
  except:
    with_background = False

  if with_background:
    try:
      backgrounddir = ast.literal_eval(cp.get('parameter_estimation', 'background'))
    except:
      with_background = False

  # check whether to plot all joint posteriors
  try:
    allposteriors = cp.getboolean('plotting', 'all_posteriors')
  except:
    allposteriors = False # default to False

  figformat = ['png'] # default to outputting png versions of figures
  # check whether to (also) output figures as eps
  try:
    witheps = cp.getboolean('plotting', 'eps_output')
  except:
    witheps = False

  if witheps:
    figformat.append('eps')

  # check whether to (also) output figures as pdf
  try:
    withpdf = cp.getboolean('plotting', 'pdf_output')
  except:
    withpdf = False

  if withpdf:
    figformat.append('pdf')

  # get the detectors to use
  ifos = [] # list of detectors to use
  withjoint = False
  preprocdat = None
  if not jointonly:
    try:
      ifos = ast.literal_eval(cp.get('general', 'detectors'))
    except:
      print("Error... could not parse list of 'detectors' in the [general] section.", file=sys.stderr)
      sys.exit(1)

    if not isinstance(ifos, list):
      print("Error... the 'detectors' value in the [general] section must be a list.", file=sys.stderr)
      sys.exit(1)

    if len(ifos) < 1:
      print("Error... the 'detectors' value in the [general] section must contain at least one detector name.", file=sys.stderr)
      sys.exit(1)

    # check whether to include the joint (multi-detector) analysis
    if len(ifos) > 1: # only required if there's been more than one detector
      try:
        withjoint = cp.getboolean('general', 'with_joint')
      except:
        withjoint = True # default to including the joint (multi-detector) analysis results

    # get paths to pre-processed data files
    try:
      preprocdat = ast.literal_eval(cp.get('data', 'files'))
    except:
      print("Error... could not parse dictionary of 'files' in the [data] section.", file=sys.stderr)
      sys.exit(1)

    if not isinstance(preprocdat, dict):
      print("Error... the 'files' value in the [data] section must be a dictionary.", file=sys.stderr)
      sys.exit(1)

    # check there is a value for each detector
    for ifo in ifos:
      if ifo not in preprocdat: # check if detector is in dictionary
        print("Error... no pre-processed data file is given for '%s'." % ifo, file=sys.stderr)
        sys.exit(1)

    # make table with all stats and plots
    datatable = create_data_table(preprocdat, outdir, figformats=figformat, harmonics=harmonics)

  # add 'Joint' on to the list of detectors if required
  if withjoint or jointonly:
    ifos.append('Joint')

  # get the posterior sample files
  try:
    postfiles = ast.literal_eval(cp.get('parameter_estimation', 'posteriors'))
  except:
    print("Error... could not parse dictionary of 'posteriors' in the [parameter_estimation] section.", file=sys.stderr)
    sys.exit(1)

  if not isinstance(postfiles, dict):
    print("Error... the 'posteriors' value in the [parameter_estimation] section must be a dictionary.", file=sys.stderr)
    sys.exit(1)

  for ifo in ifos:
    if ifo not in postfiles: # check if detector is in dictionary
      print("Error... no posterior file is given for '%s'." % ifo, file=sys.stderr)
      sys.exit(1)
    else: # check file exists
      if not os.path.isfile(postfiles[ifo]):
        print("Error... posterior file '%s' for '%s' does not exist." % (postfiles[ifo], ifo), file=sys.stderr)
        sys.exit(1)

  # create list of links to the various parts of the page
  linktable = htmltag('div', tagstyle="text-align: left; float: left")
  linkstable = htmltable()
  linkstable.addrow()

  # dictionary for output html page
  htmlinput = {}

  htmlinput['psrname'] = pname
  htmlinput['title'] = pname
  if injection:
    titlename = 'INJ ' + pname
  else:
    titlename = 'PSR ' + pname
  if atnfurl != None:
    htmlinput['h1title'] = atag(atnfurl, linktext=titlename).text
  else:
    htmlinput['h1title'] = titlename

  # create table of pulsar info from par file
  psrtable = create_psr_table(par)
  htmlinput['pulsartable'] = psrtable

  # get posterior class (containing samples, sample plots and posterior plots)
  postinfo = posteriors(postfiles, outdir, harmonics=harmonics, modeltype=modeltype,
                        biaxial=biaxial, injection=injectionfile, usegwphase=usegwphase)

  # create table of upper limits, SNR and evidence ratios
  htmlinput['limitstable'] = postinfo.create_limits_table(f0 , sdlim=sdlim, dist=dist, ul=upperlim)

  htmlinput['selectedposteriors'] = postinfo.create_joint_plots_table(title='Selected parameters')

  if allposteriors:
    htmlinput['jointposteriors'] = postinfo.create_joint_plots_table(allparams=True)
    linkstable.adddata("Posteriors ("+ atag('#selectedposteriors', linktext="Selected").text+", "+atag('#jointposteriors', linktext="All").text+")", dataclass="rightborder")
  else:
    htmlinput['jointposteriors'] = ""
    linkstable.adddata(atag('#selectedposteriors', linktext="Posteriors").text, dataclass="rightborder")

  if with_background:
    bginfo = create_background(backgrounddir, postinfo.snrs, postinfo.bsn, outdir, Bci=postinfo.bci, Bcin=postinfo.bcin)

    # add background evidence plots
    htmlinput['evidenceplots'] = bginfo.create_background_table()
    linkstable.adddata(atag("#evidenceplots", linktext="Backgrounds").text, dataclass="rightborder")
  else:
    htmlinput['evidenceplots'] = ""

  # add data plots
  htmlinput['dataplots'] = str(datatable)
  linkstable.adddata(atag("#dataplots", linktext="Data Plots").text, dataclass="rightborder")

  # add posterior samples plots
  htmlinput['posteriorsamples'] = postinfo.create_sample_plot_table(figformats=figformat)
  linkstable.adddata(atag("#posteriorsamples", linktext="Samples").text, dataclass="rightborder")

  # add statistics
  htmlinput['posteriorstats'] = postinfo.create_stats_table(credints=credints)
  linkstable.adddata(atag("#posteriorstats", linktext="Statistics").text)

  linktable.set_tagtext(linkstable.tabletext)

  # set index page link
  if indexpage != None:
    indexlink = htmltag('div', tagstyle="text-align: left; float: right; padding-right: 8px")
    indexlink.set_tagtext(atag(indexpage, linktext="Index Page").text)
    indexlinktxt = indexlink.text
  else:
    indexlinktxt = ""

  htmlinput['linkstable'] = linktable.text + indexlinktxt

  # output CSS file
  cssfile = os.path.join(outdir, 'resultspage.css')
  fp = open(cssfile, 'w')
  fp.write(result_page_css)
  fp.close()

  htmlinput['cssfile'] = os.path.basename(cssfile)

  # get time/date for file creation
  now = datetime.datetime.now()

  # add footer containing author, date and command lines used for page
  htmlinput['footer'] = "{} - {}<br><br>Command lines used:<br>{}<br>{}<br>".format(__author__, now.strftime('%a %d %b %Y'), ' '.join(sys.argv), __version__)

  # create page
  try:
    htmlfile = os.path.join(outdir, pname+'.html')
    fp = open(htmlfile, 'w')
    fp.write(htmlpage.format(**htmlinput))
    fp.close()
  except:
    print("Error... there was a problem outputting the html page.", file=sys.stderr)
    sys.exit(1)

  # output JSON file with results information
  info = {} # a dictionary with the analysis information

  info['PSR'] = pname # set the pulsar name

  # data about the pulsar
  info['Pulsar data'] = {}
  info['Pulsar data']['F0'] = f0
  info['Pulsar data']['F0ROT'] = f0
  info['Pulsar data']['F0GW'] = 2.*f0 # assuming l=m=2 mode emission
  info['Pulsar data']['F1'] = f1
  info['Pulsar data']['F1ROT'] = f1
  info['Pulsar data']['F1GW'] = 2.*f1 # assuming l=m=2 mode emission
  info['Pulsar data']['F1SD'] = f1sd # acceleration corrected f1
  info['Pulsar data']['DIST'] = dist
  info['Pulsar data']['RA'] = par['RA_RAD']   # right ascension in radians
  info['Pulsar data']['DEC'] = par['DEC_RAD'] # declination in radians
  info['Pulsar data']['START'] = par['START']
  info['Pulsar data']['FINISH'] = par['FINISH']
  info['Pulsar data']['BINARY'] = par['BINARY']
  info['Pulsar data']['PEPOCH'] = par['PEPOCH']
  info['Pulsar data']['POSEPOCH'] = par['POSEPOCH']
  info['Pulsar data']['EPHEM'] = par['EPHEM']
  info['Pulsar data']['UNITS'] = par['UNITS']
  info['Pulsar data']['ASSOC'] = assoc # any association for the pulsar
  info['Pulsar data']['spin-down limit'] = sdlim
  info['Pulsar data']['par file'] = parfile
  info['Pulsar data']['ATNF URL'] = atnfurl

  # data about the upper limits
  for ifo in ifos:
    info[ifo] = {}
    info[ifo]['Upper limits'] = {}
    info[ifo]['Upper limits']['credible region'] = upperlim
    info[ifo]['Upper limits']['H0'] = postinfo.h0_ul(ifo)
    info[ifo]['Upper limits']['ELL'] = postinfo.ellipticity_ul(ifo)
    info[ifo]['Upper limits']['Q22'] = postinfo.q22_ul(ifo)
    info[ifo]['Upper limits']['C21'] = postinfo.C21_ul(ifo)
    info[ifo]['Upper limits']['C22'] = postinfo.C22_ul(ifo)
    info[ifo]['Upper limits']['I21'] = postinfo.I21_ul(ifo)
    info[ifo]['Upper limits']['I31'] = postinfo.I21_ul(ifo)
    info[ifo]['Upper limits']['spin-down ratio'] = postinfo.sdlim_ratio(ifo)

    # evidence ratios
    info[ifo]['Bayes factors'] = {}
    info[ifo]['Bayes factors']['Signal vs Noise'] = postinfo.bsn[ifo]

    if 'Joint' == ifo:
      info[ifo]['Bayes factors']['Coherent vs Incoherent'] = postinfo.bci
      info[ifo]['Bayes factors']['Coherent vs Incoherent or Noise'] = postinfo.bcin

    info[ifo]['SNR'] = postinfo.snr(ifo)

    # amplitude noise spectal densities
    if ifo != 'Joint':
      info[ifo]['Amplitude spectral density'] = {}

      if len(harmonics) == 1:
        info[ifo]['Amplitude spectral density']['Spectrum'] = datatable.asds[harmonics[0]][ifo].tolist()
        info[ifo]['Amplitude spectral density']['Mean'] = np.mean(datatable.asds[harmonics[0]][ifo])
        info[ifo]['Amplitude spectral density']['Median'] = np.median(datatable.asds[harmonics[0]][ifo])
        info[ifo]['Amplitude spectral density']['Maximum'] = np.max(datatable.asds[harmonics[0]][ifo])
      else:
        for h in harmonics:
          info[ifo]['Amplitude spectral density']['%df harmonic' % int(h)] = {}
          info[ifo]['Amplitude spectral density']['%df harmonic' % int(h)]['Spectrum'] = datatable.asds[h][ifo].tolist()
          info[ifo]['Amplitude spectral density']['%df harmonic' % int(h)]['Mean'] = np.mean(datatable.asds[h][ifo])
          info[ifo]['Amplitude spectral density']['%df harmonic' % int(h)]['Median'] = np.median(datatable.asds[h][ifo])
          info[ifo]['Amplitude spectral density']['%df harmonic' % int(h)]['Maximum'] = np.max(datatable.asds[h][ifo])

  jsonfile = os.path.join(outdir, pname+'.json')
  fp = open(jsonfile, 'w')
  json.dump(info, fp, indent=2)
  fp.close()

  sys.exit(0)
