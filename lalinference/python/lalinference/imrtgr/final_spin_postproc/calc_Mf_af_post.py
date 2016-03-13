# Code to read in LALInference posterior samples and plot the histogram of the posteriors on the final mass and spin (in addition to outputting medians and confidence intervals), using various fits, some of which include the effects of in-plane spins on the final spin
# NKJ-M, 03.2016, based on earlier code

import numpy as np
import argparse

from Mfaf_posterior import *

# Set up the parsing

parser = argparse.ArgumentParser(description = 'Calculate the posteriors on the radiated mass and peak luminosity from LALInference posterior samples using a variety of NR fits for these quantities.')
parser.add_argument("datfile", help = "posterior_samples.dat file to read in")
parser.add_argument("tag", help = "tag for output files")
parser.add_argument("outdir", help = "directory for output file (default: current directory)", default = ".")
parser.add_argument("--outputsamples", help = "output posterior samples (currently just for the HLZ Erad and the SHKJDK6mode Lpeak)", action = "store_true")
args = parser.parse_args()

print(args.tag)
print("%s\n"%(args.outdir))

# read the posterior samples for the masses and spins
data = np.genfromtxt(args.datfile, dtype=None, names=True)
if ('m1_source' in data.dtype.names) and ('m2_source' in data.dtype.names):
  m1, m2 = data['m1_source'], data['m2_source']
else:
  print("WARNING! Source frame masses not available, using redshifted masses instead.")
  m1, m2 = data['m1'], data['m2']
if ('a1' in data.dtype.names) and ('a2' in data.dtype.names):
  chi1, chi2 = data['a1'], data['a2']
else:
  print("WARNING! Spins not found--setting to zero.")
  chi1, chi2 = np.zeros(len(m1)), np.zeros(len(m2))
if ('tilt1' in data.dtype.names) and ('tilt2' in data.dtype.names):
  tilt1, tilt2 = data['tilt1'], data['tilt2']
else:
  print("WARNING! tilts not found--using aligned spin components.")
  chi1, chi2 = data['a1z'], data['a2z']
  tilt1, tilt2 = np.zeros(len(m1)), np.zeros(len(m2))
if 'chi_p' in data.dtype.names:
  chi_p = data['chi_p']
else:
  chi_p = np.zeros(len(m1))
if 'phi12' in data.dtype.names:
  phi12 = data['phi12']
else:
  phi12 = np.zeros(len(m1))

# Loop over the final mass (and spin) fits

# Select the appropriate fits (the IMRPhenomPv2 fit requires chi_p, which we do not have for the evolved spin samples)

if np.count_nonzero(chi_p) > 0:
  fits = ['HLZ', 'HLZinplane', 'HLZinplaneIMRPhenomP', 'Husaetal', 'IMRPhenomPv2', 'IMRPhenomPtilts', 'BRaligned', 'BR']
else:
  fits = ['HLZ', 'HLZinplane', 'HLZinplaneIMRPhenomP', 'Husaetal', 'IMRPhenomPtilts', 'BRaligned', 'BR']

for fitsel in fits:
  print("Calculating for fit %s"%(fitsel))
  if fitsel == 'HLZinplaneIMRPhenomP':
    Mf_samples, af_samples = Mfaf_posterior(fitsel, m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12, args.tag, args.outdir)
  else:
    Mfaf_posterior(fitsel, m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12, args.tag, args.outdir)

if args.outputsamples:
  np.savetxt("%s/%s_Mfaf_samples.dat"%(args.outdir, args.tag), zip(Mf_samples, af_samples), header = "Mf_source\t af", delimiter = "\t")
