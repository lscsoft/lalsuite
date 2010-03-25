#!/usr/bin/env python

import nds2
import pylab
import numpy

daq = nds2.Daq('ldas-pcdev1.ligo.caltech.edu',31200)

gps_start = 940523050
gps_end   = 940523060

# Get three channels in one request...
first_results = daq.get_data(
     ['L0:PEM-EX_SEISX', 'L0:PEM-EX_SEISY', 'L0:PEM-EX_SEISZ'], 
     gps_start, 
     gps_end)

pylab.figure(0)

for result in first_results:
    dt    = 1.0 / result.signal_rate
    times = numpy.arange(gps_start, gps_end, dt)

    pylab.plot(times, result.data, label=result.name)

pylab.legend()
pylab.title("EX Seismic")



# Get three more channels... without needing to shut down first!
second_results = daq.get_data(
    ['L0:PEM-EY_SEISX', 'L0:PEM-EY_SEISY', 'L0:PEM-EY_SEISZ'], 
    gps_start, 
    gps_end)

pylab.figure(1)

for result in second_results:
    dt    = 1.0 / result.signal_rate
    times = numpy.arange(gps_start, gps_end, dt)

    pylab.plot(times, result.data, label=result.name)

pylab.legend()
pylab.title("EY Seismic")

pylab.show()


