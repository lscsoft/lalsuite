# Code to read in LALInference posterior samples and plot the histogram of the posteriors on the final spin (in addition to outputting medians and confidence intervals), using various fits, some of which include the effects of in-plane spins on the final spin
# NKJ-M, 03.2016, based on earlier code

import numpy as np
import argparse
import sys

from af_posterior import *

# Set up the parsing

parser = argparse.ArgumentParser(description = 'Calculate the posteriors on the radiated mass and peak luminosity from LALInference posterior samples using a variety of NR fits for these quantities.')
parser.add_argument("datfile", help = "posterior_samples.dat file to read in")
parser.add_argument("tag", help = "tag for output files")
parser.add_argument("outdir", help = "directory for output file (default: current directory)", default = ".")
parser.add_argument("--outputsamples", help = "output posterior samples (currently just outputs the best fit)", action = "store_true")
args = parser.parse_args()

print(args.tag)
print("%s\n"%(args.outdir))

# read the posterior samples for the masses and spins (we can just use the detector frame masses here, since we're just computing the final spin, which only needs the mass ratio)
data = np.genfromtxt(args.datfile, dtype=None, names=True)
if ('m1_source' in data.dtype.names) and ('m2_source' in data.dtype.names):
  m1, m2 = data['m1_source'], data['m2_source']
elif ('m1' in data.dtype.names) and ('m2' in data.dtype.names):
  m1, m2 = data['m1'], data['m2']
else:
  print("FAILURE! Masses not found--there is likely something wrong with the posterior file you input.")
  sys.exit()
if not (('a1' in data.dtype.names) and ('a2' in data.dtype.names)):
  print("WARNING! Spins not found--setting to zero.")
  chi1, chi2 = np.zeros(len(m1)), np.zeros(len(m2))
elif ('tilt1' in data.dtype.names) and ('tilt2' in data.dtype.names):
  chi1, chi2 = data['a1'], data['a2']
  tilt1, tilt2 = data['tilt1'], data['tilt2']
else:
  print("WARNING! tilts not found--using aligned spin components.")
  chi1, chi2 = data['a1z'], data['a2z']
  tilt1, tilt2 = np.zeros(len(m1)), np.zeros(len(m2))
if 'chi_p' in data.dtype.names:
  chi_p = data['chi_p']
else: # Compute chi_p using the expression in Hannam et al. PRL 113, 151101 (2014) [see also Eq. (29) in https://dcc.ligo.org/LIGO-T1500602]; note that those expressions use the opposite mass ordering convention to LALInference (i.e., m2 \geq m1 instead of m1 \geq m2)                                                
  print("NOTE: chi_p not found; Computing chi_p from spins and tilts.")
  A1 = 2. + 1.5*m2/m1
  A2 = 2. + 1.5*m1/m2
  chi_p = np.maximum(chi1*np.sin(tilt1),(A2/A1)*(m2*m2/(m1*m1))*chi2*np.sin(tilt2))
if 'phi12' in data.dtype.names:
  phi12 = data['phi12']
else:
  print("WARNING! phi12 not found--setting to zero.")
  phi12 = np.zeros(len(m1))

# Loop over the final spin fits

fits = ['HLZ', 'HLZinplane', 'HLZinplaneIMRPhenomP', 'Husaetal', 'IMRPhenomPv2', 'IMRPhenomPtilts', 'BRaligned', 'BR']

for fitsel in fits:
  print("Calculating for fit %s"%(fitsel))
  if fitsel == 'HLZinplaneIMRPhenomP':
    af_samples = af_posterior(fitsel, m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12, args.tag, args.outdir)
  else:
    af_posterior(fitsel, m1, m2, chi1, chi2, tilt1, tilt2, chi_p, phi12, args.tag, args.outdir)

if args.outputsamples:
  np.savetxt("%s/%s_af_samples.dat"%(args.outdir, args.tag), af_samples, header = "af")
