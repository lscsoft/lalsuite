#!/usr/bin/env python
"""

Demo of pyNDS, showing online seismic channels from LHO.

"""
__author__       = "Leo Singer <leo.singer@ligo.org>"
__organization__ = ["LIGO", "California Institute of Technology"]
__copyright__    = "Copyright 2010, Leo Singer"



import nds
from pylab import *


# Connect to the NDS1 server at Hanford.
# Note that since no protocol is specified, it will first attempt NDS2, then
# fall back to NDS1.
daq = nds.daq('blue.ligo-wa.caltech.edu', 31200)


# Recieve list of 'online' channels
channels = daq.recv_channel_list(nds.channel_type.online)


# Get the seismic channels at EX
desired_channels = ('H0:PEM-EY_SEISX', 'H0:PEM-EY_SEISY', 'H0:PEM-EY_SEISZ')
daq.request_channels(c for c in channels if c.name in desired_channels)


# Create a new figure
fig = figure()
fig.show()


# Iterate over blocks of data as they com in
for data in daq.seek(0, 0, 1):
    
    # Get timestamp (tuple of GPS seconds and nanoseconds)
    timestamp = daq.timestamp[0] + daq.timestamp[1]*0.000000001
    
    # Clear axes
    cla()
    for i in range(3):
        # Compute time scale
        t = arange(len(data[i])) * 1000.0 / daq.requested_channels[i].rate
        # Plot
        plot(t, data[i], label=str(daq.requested_channels[i]))
    
    ylim((-6000, 6000))
    xlabel('milliseconds relative to GPS time %.1f' % timestamp)
    legend()
    draw()
