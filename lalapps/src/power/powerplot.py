#!/usr/bin/env python
import sys
from matplotlib.patches import Patch
from pylab import *

class C: pass

def get_burst_table(triggerfile):
    M = load('%s' % triggerfile)
    c = C()
    c.start_time = M[:,0]+1.0e-9*M[:,1]
    c.peak_time = M[:,2]+1.0e-9*M[:,3]
    c.duration = M[:,4]
    c.central_freq = M[:,5]
    c.bandwidth = M[:,6]
    c.snr = M[:,7]
    c.confidence = M[:,8]
    return c

#e12 = get_burst_table('data/test.txt')
