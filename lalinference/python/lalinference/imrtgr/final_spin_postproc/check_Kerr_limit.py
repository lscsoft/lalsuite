# Check that the IMRPhenomP precessing extension to various aligned-spin final spin fits obeys the Kerr limit using a mesh. We just check this for various mass ratios and tilt angles for extremal initial spins, since the maximum in-plane spin contribution for given aligned spin components comes from taking the in-plane spins to have the maximum value and be aligned with each other.
# NKJ-M, 03.2016

import numpy as np
import nrutils_new as nr

# Set up mass ratio range to consider and mesh

qmin = 1. #1. # Minimum mass ratio (q \geq 1)
qmax = 50. #20 #100. # Maximum mass ratio (q \geq 1)
nq = 83 #20 #137 #51 #27 #100 # Number of q points
ntilt = 107 #20 #137 #53 #106 #73 #217 #100 # Number of tilt points

# Do we plot? (Only for the mesh version...)

plot = 1

# Set up variables and compute

chifmax = 0.

chifqmax = np.zeros(nq)
chifqmaxtilt1 = np.zeros(nq)
chifqmaxtilt2 = np.zeros(nq)

i = 0

qvec = np.linspace(qmin,qmax,num=nq,endpoint=True)

for q in qvec:
    print("q = %f"%(q))
    for tilt1 in np.linspace(0.,np.pi,num=ntilt,endpoint=True):
        for tilt2 in np.linspace(0.,np.pi,num=ntilt,endpoint=True):
            chif = nr.bbh_final_spin_precessing_Healyetal_extension_Minit(q,1.,1.,1.,tilt1,tilt2,0.)
            #chif = nr.bbh_final_spin_precessing_Healyetal_extension_Mf(q,1.,1.,1.,tilt1,tilt2,0.)
            #chif = nr.bbh_final_spin_precessing_Husaetal_extension_Minit(q,1.,1.,1.,tilt1,tilt2,0.)
            if abs(chif) > chifqmax[i]:
                chifqmax[i] = abs(chif)
                chifqmaxtilt1[i], chifqmaxtilt2[i] = tilt1, tilt2
                if chifqmax[i] > chifmax:
                    chifmax = chifqmax[i]
                    chifmaxq, chifmaxtilt1, chifmaxtilt2 = q, tilt1, tilt2
    print("Maximum chif for q = %f with tilt1 = %f and tilt2 = %f: %f"%(q, chifqmaxtilt1[i], chifqmaxtilt2[i], chifqmax[i]))
    i += 1

print("Overall maximum chif of %f for q = %f, tilt1 = %f, tilt2 = %f"%(chifmax, chifmaxq, chifmaxtilt1, chifmaxtilt2))

if plot:
    import matplotlib.pyplot as plt
    plt.plot(qvec, chifqmax)
    plt.xlabel('$q$')
    plt.ylabel('$\chi_f^\mathrm{max}$')
    plt.show()
    plt.close()
