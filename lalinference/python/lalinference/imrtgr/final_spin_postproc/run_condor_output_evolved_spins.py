'''
Code to read in LALInference posterior samples at each node, evolve the spins up to the ISCO, and output the resulting mass and spin samples. This code takes in the input posterior sample file with name `sometag.nodenum` and outputs `..._nodenum.dat`.

Usage:
python run_condor_output_evolved_spins.py <tag for the sub file> <output directory> <Node number>
'''
# NKJ-M, 02.2016; Anuradha Gupta, 03.2016

import numpy as np
import argparse
import pneqns
import time

# Settings
v_final = 6.**-0.5  # some fiducial value for final orbital velocity. Currently, it is Schwarzschild ISCO velocity.
flow = 20 # Start frequency in Hz
dt = 1./5120. # steps in time for the integration

# Set up the parsing
parser = argparse.ArgumentParser(description = 'Evolve the LALInference posterior samples of spins')
parser.add_argument("tag", help = "tag for output files")
parser.add_argument("outdir", help = "directory for output file (default: current directory)", default = ".")
parser.add_argument("nodenum", help = "Node number for the condor runs", type=int)
args = parser.parse_args()

infile = args.tag+'.%s'%(args.nodenum)

# read the posterior samples for the masses and spins
data = np.genfromtxt(infile, dtype=None, names=True)
if ('m1_source' in data.dtype.names) and ('m2_source' in data.dtype.names):
  m1, m2 = data['m1_source'], data['m2_source']
else:
  print("WARNING! Source frame masses not available, using redshifted masses instead.")
  m1, m2 = data['m1'], data['m2']
if ('a1' in data.dtype.names) and ('a2' in data.dtype.names):
  chi1, chi2 = data['a1'], data['a2']
else:
  chi1, chi2 = np.zeros(len(m1)), np.zeros(len(m2))
if ('tilt1' in data.dtype.names) and ('tilt2' in data.dtype.names):
  tilt1, tilt2 = data['tilt1'], data['tilt2']
else:
  chi1, chi2 = data['a1z'], data['a2z']
  tilt1, tilt2 = np.zeros(len(m1)), np.zeros(len(m2))
if 'phi12' in data.dtype.names:
  phi12 = data['phi12']
else:
  phi12 = np.zeros(len(m1))

# Evolve spins
Msun2s = 4.92549e-6 # s/Msun


start = time.clock()

for i in range(len(m1)):
    v0 = ((m1[i]+m2[i])*Msun2s*np.pi*flow)**(1./3)
    print "Initial values of angles:"
    print "cos tilt1 = ", np.cos(tilt1[i]), "cos tilt2 = ", np.cos(tilt2[i]), "cos phi12 = ", np.cos(phi12[i])
    if np.logical_and(tilt1[i] < 1e-3, tilt2[i] < 1e-3):
        phi12[i] = 0.
    else:
        tilt1[i], tilt2[i], phi12[i] = pneqns.find_tilts_and_phi12_at_freq(v0, m1[i], m2[i], chi1[i]*np.sin(tilt1[i]), 0., chi1[i]*np.cos(tilt1[i]), chi2[i]*np.sin(tilt2[i])*np.cos(phi12[i]), chi2[i]*np.sin(tilt2[i])*np.sin(phi12[i]), chi2[i]*np.cos(tilt2[i]), 0., 0., 1., v_final, dt)

print("Elapsed time for %f calculations: %f s"%(i, time.clock()-start))

np.savetxt("%s/%s_evolved_spin_samples_%s.dat"%(args.outdir, args.tag, args.nodenum), zip(m1, m2, chi1, chi2, tilt1, tilt2, phi12), header = "m1_source\t m2_source\t a1\t a2\t tilt1\t tilt2\t phi12", delimiter = "\t")
