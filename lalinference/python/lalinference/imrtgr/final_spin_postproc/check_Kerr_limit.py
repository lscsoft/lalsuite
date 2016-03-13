# Check that the IMRPhenomP precessing extension to various aligned-spin final spin fits obeys the Kerr limit. We just check this for various mass ratios and tilt angles for extremal initial spins, since the maximum in-plane spin contribution for given aligned spin components comes from taking the in-plane spins to have the maximum value and be aligned with each other.
# NKJ-M, 03.2016

import numpy as np
import nrutils_new as nr

# Set up mass ratio range to consider and mesh

qmin = 1. #1. # Minimum mass ratio (q \geq 1)
qmax = 20. #100. # Maximum mass ratio (q \geq 1)
nq = 27 #100 # Number of q points
ntilt = 53 #106 #73 #217 #100 # Number of tilt points

chifmax = 0.

for q in np.linspace(qmin,qmax,num=nq,endpoint=True):
    print("q = %f"%(q))
    chifqmax = 0.
    for tilt1 in np.linspace(0.,np.pi,num=ntilt,endpoint=True):
        for tilt2 in np.linspace(0.,np.pi,num=ntilt,endpoint=True):
            chif = nr.bbh_final_spin_precessing_Healyetal_extension_Minit(q,1.,1.,1.,tilt1,tilt2,0.)
            #chif = nr.bbh_final_spin_precessing_Healyetal_extension_Mf(q,1.,1.,1.,tilt1,tilt2,0.)
            #chif = nr.bbh_final_spin_precessing_Husaetal_extension_Minit(q,1.,1.,1.,tilt1,tilt2,0.)
            if abs(chif) > chifqmax:
                chifqmax = abs(chif)
                chifqmaxtilt1, chifqmaxtilt2 = tilt1, tilt2
                if chifqmax > chifmax:
                    chifmax = chifqmax
                    chifmaxq, chifmaxtilt1, chifmaxtilt2 = q, tilt1, tilt2
    print("Maximum chif for q = %f with tilt1 = %f and tilt2 = %f: %f"%(q, chifqmaxtilt1, chifmaxtilt2, chifqmax))
            

print("Overall maximum chif of %f for q = %f, tilt1 = %f, tilt2 = %f"%(chifmax, chifmaxq, chifmaxtilt1, chifmaxtilt2))
