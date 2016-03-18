# Check that the IMRPhenomP precessing extension to various aligned-spin final spin fits obeys the Kerr limit using a random sampling. We just check this for various mass ratios and tilt angles for extremal initial spins, since the maximum in-plane spin contribution for given aligned spin components comes from taking the in-plane spins to have the maximum value and be aligned with each other.
# NKJ-M, 03.2016

import numpy as np
import nrutils_new as nr
import time

# Set up mass ratio range to consider and number of points

qmin = 1. #1. # Minimum mass ratio (q \geq 1)
qmax = 50. #10. #100. # Maximum mass ratio (q \geq 1)

npts = 1e3

# Choose whether to sample in tilts or cos(tilts)
costilts = 0

chifmax = 0.

# Start timing

tstart = time.clock()

# Make random numbers in [0,1)

rq, rtilt1, rtilt2 = np.random.random((3,npts))

# Rescale

q = rq*(qmax-qmin) + qmin

if costilts:
    tilt1 = np.arccos(rtilt1*2.-1.)
    tilt2 = np.arccos(rtilt2*2.-1.)
else:
    tilt1 = rtilt1*np.pi
    tilt2 = rtilt2*np.pi

#chif = nr.bbh_final_spin_precessing_Healyetal_extension_Minit(q,1.,1.,1.,tilt1,tilt2,0.)
chif = nr.bbh_final_spin_precessing_Healyetal_extension_Mf(q,1.,1.,1.,tilt1,tilt2,0.)                                                                                   
#chif = nr.bbh_final_spin_precessing_Husaetal_extension_Minit(q,1.,1.,1.,tilt1,tilt2,0.) 

chifmax = np.max(chif)

idxmax = np.argmax(chif)

chifmaxq, chifmaxtilt1, chifmaxtilt2 = q[idxmax], tilt1[idxmax], tilt2[idxmax]

qmax = np.max(q)

maxqchif = chif[np.argmax(q)]

# End time

tend = time.clock()

print("Computation time for %g samples: %f s"%(npts,tend-tstart))

print("Overall maximum chif of %f for q = %f, tilt1 = %f, tilt2 = %f"%(chifmax, chifmaxq, chifmaxtilt1, chifmaxtilt2))

print("Maximum q sampled is %f, with a chif of %f"%(qmax,maxqchif))
