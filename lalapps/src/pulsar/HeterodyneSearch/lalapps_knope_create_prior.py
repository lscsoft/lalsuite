# -*- coding: utf-8 -*-
#
#       lalapps_knope_create_prior.py
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

"""
A script to create a prior file based on a given configuration file.

If the configuration (.ini) file contains parameters and ranges for them then
this will be output into the prior file.

In general the any amplitude range specified in the [parameter] section will
be used for the priors. However, if a 'priorfile' is given in the [prior] section
then an upper limit for any required amplitude parameters will be searched for
within this file. This file will be a JSON file with a dictionary keyed on pulsar
names in the following format:
{
  "J0534+2200": {"H0UL": 2.3e-25, "C22UL": None, "C21UL": None, "I31UL": None, "I21UL": None}
}
giving 95% upper limits on any possible amplitude parameters. These will be
used to place Fermi-Dirac priors on the amplitude parameters.

If the 'priorfile' is not given then the amplitude priors will be derived from
"expected" 95% upper limits from previous runs based on estimates of the amplitude
spectral densities from a "Joint" multidetector analysis. With the frequency
at which to define the priors dependent on the 'modeltype'. Again this will be
used to inform a Fermi-Dirac prior.

[pulsar]
# the pulsar name
name = J0534+2200
# the rotational frequency of the pulsar
f0 = 23.4
# a list of the frequency scale factors used in the analysis
freqfactors = [2.0]

[parameters]
# the set of parameters to use in the prior file
H0 = {"priortype": "uniform", "ranges": [0.0, 1e-22]},
PHI0 = {"priortype": "gaussian", "ranges": [1.0, 1.0]}

[prior]
# a (JSON) file containing any prior upper limits from previous runs
priorfile = S6priors.txt
# the model type 'waveform' or 'source'
modeltype = waveform
# a dictionary of paths to previous ASD files from which to derive a prior
asd_files = {'H1': path_to_S6_H1_ASD_file, 'L1': path_to_S6_H1_ASD_file, 'V1': path_to_VSR4_V1_ASD_file}
# a dictionary of observation times (days) for each detector
obs_times = {'H1': 100, 'L1': 120, 'V1': 40}
# the type of prior: 'fermidirac' or 'uniform'
priortype = fermidirac

[output]
# the output prior file
outputfile = PSRJ_H1.prior

"""

# make print statements python 3-proof
from __future__ import print_function, division

import sys
import os
import numpy as np
import json
import argparse
from ConfigParser import ConfigParser
import ast


def fermidirac_rsigma(ul, mufrac=0.4, cdf=0.95):
  """
  Calculate the r and sigma parameter of the Fermi-Dirac distribution to be used.

  Based on the definition of the distribution given in https://www.authorea.com/users/50521/articles/65214/_show_article
  the distribution will be defined by a mu parameter at which the distribution has 50% of it's maximum
  probability, and mufrac which is the fraction of mu defining the range from which the distribution falls from
  97.5% of the maximum down to 2.5%. Using an upper limit defining a given cdf of the distribution the parameters
  r and sigma will be returned.
  """

  from scipy.optimize import root

  Z = 7.33 # factor that defined the 97.5% -> 2.5% probability attenuation band around mu

  r = 0.5*Z/mufrac # set r

  # using the Fermi-Dirac CDF to find sigma given a distribution where the cdf value is found at ul
  sol = root(lambda s: cdf*np.log(1.+np.exp(r))-np.log(1.+np.exp(-r))-(ul/s)-np.log(1.+np.exp((ul/s)-r)), ul)
  sigma = sol.x[0]

  return r, sigma


# main function
if __name__=='__main__':
  description = """This script will create a prior file for use by lalapps_pulsar_parameter_estimation_nested."""

  parser = argparse.ArgumentParser( description = description )
  parser.add_argument("inifile", help="The configuration (.ini) file")

  # parse input options
  opts = parser.parse_args()

  inifile = opts.inifile

  # read in config file
  cp = ConfigParser()
  try:
    cp.optionxform = str # make parser case sensitive
    cp.read(inifile)
  except:
    print("Error... cannot parse configuration file '%s'" % inifile, file=sys.stderr)


  try:
    outfile = cp.get('output', 'outputfile')
  except:
    print("Error... no 'outputfile' specified in [output] section of the configuration file.", file=sys.stderr)
    sys.exit(1)

  # open output file for writing
  try:
    fp = open(outfile, 'w')
  except:
    print("Error... could not open output prior file '%s'" % outfile)
    sys.exit(1)

  # get parameters from [parameters] section
  parameters = {}
  if cp.has_section('parameters'):
    for item in cp.items('parameters'):
      parameters[item[0].upper()] = ast.literal_eval(item[1])
    if len(parameters) == 0:
      print("Error... no parameters given in the [parameters] section of the configuration file.", file=sys.stderr)
  else:
    print("Error... no [parameters] section in the configuration file.", file=sys.stderr)

  # Go through parameters and output to prior file
  for param in parameters:
    if 'priortype' not in parameters[param]:
      print("Error... no 'priortype' given for parameter '%s'" % param, file=sys.stderr)
      sys.exit(1)
    if 'ranges' not in parameters[param]:
      print("Error... no 'ranges' given for parameter '%s'" % param, file=sys.stderr)
      sys.exit(1)

    ptype = parameters[param]['priortype']
    ranges = parameters[param]['ranges']

    if len(ranges) != 2:
      print("Error... 'ranges' for parameter '%s' must be a list or tuple with two entries" % param, file=sys.stderr)
      sys.exit(1)

    # write to file
    fp.write('%s\t%s\t%.16le\t%.16le\n' % (param, ptype, ranges[0], ranges[1]))

  # check if a previous prior file exists from which to get amplitude prior information

  # if required use Bk files to set amplitude limits dependent on the required waveform model type
  if cp.has_section('prior'):
    try:
      modeltype = cp.get('prior', 'modeltype')
    except:
      print("Error... no 'modeltype' given in [prior] section of configuration file.", file=sys.stderr)
      sys.exit(1)

    # get the pulsar name
    try:
      pname = cp.get('pulsar', 'name')
    except:
      pname = None

    # get the pulsar frequency
    try:
      freq = cp.getfloat('pulsar', 'f0')
    except:
      freq = 0.

    # get the frequency factors
    try:
      freqfactors = ast.literal_eval(cp.get('pulsar', 'freqfactors'))
      sorted(freqfactors) # make sure they are in ascending order
    except:
      print("Error... no 'freqfactors' given in [pulsar] section of configuration file.", file=sys.stderr)
      sys.exit(1)

    # set the required amplitude priors
    requls = {} # dictionary to contain the required upper limits
    if modeltype == 'waveform':
      if 2. in freqfactors:
        requls['C22'] = None
      if 1. in freqfactors:
        requls['C21'] = None
    elif modeltype == 'source':
      if len(freqfactors) == 1:
        requls['H0'] = None
      if len(freqfactors) == 2:
        if 1. in freqfactors and 2. in freqfactors:
          requls['I21'] = None
          requls['I31'] = None

    if len(requls) == 0:
      print("Error... unknown frequency factors or model type in configuration file.", file=sys.stderr)
      sys.exit(1)

    # try and get a prior file
    try:
      priorfile = cp.get('prior', 'priorfile')
    except:
      priorfile = ""
    if len(priorfile) > 1:
      # check file can be read
      try:
        fpp = open(priorfile, 'r')
        priorinfo = json.load(fpp) # should be JSON file
        fpp.close()
      except:
        print("Error... could not parse prior file '%s'." % priorfile, file=sys.stderr)
        sys.exit(1)

      # see if pulsar is in prior file
      if pname in priorinfo:
        uls = priorinfo[pname]
        for ult in requls:
          if ult == 'C22':
            if 'C22UL' not in uls and 'H0UL' in uls:
              # use 'H0' value for 'C22' if present
              requls['C22'] = uls['H0UL']
          else:
            if ult+'UL' in uls:
              requls[ult] = uls[ult+'UL']

    # if there are some required amplitude limits that have not been obtained try and get amplitude spectral densities
    if None in requls.values() and freq > 0.0:
      try:
        try:
          asdfiles = ast.literal_eval(cp.get('prior', 'asd_files'))
        except:
          asdfiles = cp.get('prior', 'asd_files')
      except:
        print("Error... no ASD files given in configuration file.", file=sys.stderr)
        sys.exit(1)

      # get observation times
      try:
        try:
          obstimes = ast.literal_eval(cp.get('prior', 'obs_times'))
        except:
          obstimes = cp.get('prior', 'obs_times')
      except:
        print("Error... no observation times given in configuration file.", file=sys.stderr)
        sys.exit(1)

      if not isinstance(asdfiles, dict): # if a single file is given convert into dictionary
        asdfilestmp = {}
        obstimestmp = {}
        asdfilestmp['det'] = asdfiles
        obstimestmp['det'] = float(obstimes)
        asdfiles = asdfilestmp
        obstimes = obstimestmp

      freqmax = np.inf # get the 'minimum' maximum frequency
      freqmin = 0.     # get the 'maximum' minimum frequency
      asdlist = []
      if isinstance(asdfiles, dict): # read in all files
        for dk in asdfiles:
          if dk not in obstimes:
            print("Error... no corresponding observation times for detector '%s'" % dk, file=sys.stderr)
            sys.exit(1)
          else:
            if not isinstance(obstimes[dk], float) and not isinstance(obstimes[dk], int):
              print("Error... observation time must be a float or int.", file=sys.stderr)
              sys.exit(1)
            if not os.path.isfile(asdfiles[dk]):
              print("Error... ASD file '%s' does not exist." % asdfiles[dk], file=sys.stderr)
              sys.exit(1)
            else:
              try:
                asd = np.loadtxt(asdfiles[dk], comments=['%', '#'])
                # set min and max frequencies
                asdv = [] # empty array
                if 1. in freqfactors and (asd[0,0] <= freq and asd[-1,0] >= freq): # add ASD at 1f
                  idxf = (np.abs(asd[:,0]-freq)).argmin() # get value nearest required frequency
                  asdv.append(asd[idxf,1])
                if 2. in freqfactors and (asd[0,0] <= 2.*freq and asd[-1,0] >= 2.*freq):
                  idxf = (np.abs(asd[:,0]-2.*freq)).argmin() # get value nearest required frequency
                  asdv.append(asd[idxf,1])

                if len(asdv) > 0:
                  asdlist.append(np.array(asdv)**2/(obstimes[dk]*86400.))
                else:
                  print("Error... frequency range in ASD file does not span pulsar frequency.", file=sys.stderr)
                  sys.exit(1)
              except:
                print("Error... could not load file '%s'." % asdfiles[dk], file=sys.stderr)
                sys.exit(1)

      # get upper limit spectrum (harmonic mean of all the weighted spectra)
      mspec = np.zeros(len(freqfactors))
      for asdv in asdlist:
        # interpolate frequencies
        mspec = mspec + (1./asdv)

      mspec = np.sqrt(1./mspec) # final weighted spectrum
      ulspec = 10.8*mspec # scaled to given "averaged" 95% upper limit estimate

      # set upper limits for creating priors
      if modeltype == 'waveform':
        if 1. in freqfactors:
          if requls['C21'] == None:
            requls['C21'] = ulspec[freqfactors.index(1.0)]
        if 2. in freqfactors:
          if requls['C22'] == None:
            requls['C22'] = ulspec[freqfactors.index(2.0)]
      if modeltype == 'source':
        if len(freqfactors) == 1:
          if requls['H0'] == None:
            requls['H0'] = ulspec[0]
        else:
          if 1. in freqfactors and 2. in freqfactors:
            # set both I21 and I31 to use the maximum of the 1f and 2f es
            if requls['I21'] == None:
              requls['I21'] = np.max(ulspec)
            if requls['I31'] == None:
              requls['I31'] = np.max(ulspec)

    # get prior type
    try:
      priortype = cp.get('prior', 'priortype')
    except:
      print("Error... must specify 'priortype' in [prior] section.", file=sys.stderr)
      sys.exit(1)

    if priortype not in ['fermidirac', 'uniform']:
      print("Error... 'priortype' must be 'fermidirac' or 'uniform'", file=sys.stderr)
      sys.exit(1)

    # go through required upper limits and output a Fermi-Dirac prior that also has a 95% limit at that value
    for ult in requls:
      if requls[ult] == None:
        print("Error... a required upper limit for '%s' is not available." % ult, file=sys.stderr)
        sys.exit(1)
      else:
        if priortype == 'fermidirac':
          try:
            b, a = fermidirac_rsigma(requls[ult])
          except:
            print("Error... problem deriving the Fermi-Dirac prior for '%s'." % ult, file=sys.stderr)
            sys.exit(1)
        else:
          a = 0. # uniform prior bound at 0
          b = requls[utl]/0.95 # stretch limit to ~100% bound

        fp.write('%s\t%s\t%.16le\t%.16le\n' % (ult, priortype, a, b))

  fp.close()
  sys.exit(0)
