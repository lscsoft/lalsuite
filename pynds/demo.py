#!/usr/bin/env python
"""

Demo of pyNDS, showing archived seismic channels from LHO

"""
__author__       = "Leo Singer <leo.singer@ligo.org>"
__organization__ = ["LIGO", "California Institute of Technology"]
__copyright__    = "Copyright 2010, Leo Singer"



import nds
from pylab import *


# Connect to the NDS2 server at Hanford.
daq = nds.daq('ldas-pcdev1.ligo.caltech.edu', 31200)


# Recieve list of 'raw' channels
channels = daq.recv_channel_list(nds.channel_type.raw)


# Get the seismic channels at EX
desired_channels = ('H0:PEM-EY_SEISX', 'H0:PEM-EY_SEISY', 'H0:PEM-EY_SEISZ')
daq.request_channels(c for c in channels if c.name in desired_channels)


# Create new figure
fig = figure()
fig.show()


# Iterate over blocks of data as they come in
start = 953618739
for data in daq.seek(start, start + 3600, 60):
    
    # Get timestamp (tuple of GPS seconds and nanoseconds)
    timestamp = daq.timestamp[0] + daq.timestamp[1]*0.000000001
    
    cla()
    for i in range(3):
        # Compute time scale
        t = arange(len(data[i])) / daq.requested_channels[i].rate
        # Plot
        plot(t, data[i], label=str(daq.requested_channels[i]))
    
    xlabel('seconds relative to GPS time %.1f' % timestamp)
    ylim((-6000, 6000))
    legend()
    draw()
