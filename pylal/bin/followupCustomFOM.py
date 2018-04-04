#!/usr/bin/env python
#
# Copyright (C) 2009 Cristina Valeria Torres
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#
"""
This script is intended to recreate some basic
Figure Of Merit plots for following up specific
GPS times with cleaner plots instead of the
RoboScimon plots in the elog.
"""
__author__  = "Cristina Valeria Torres <cristina.torres@ligo.org>"
__prog__    = 'followupCustomFOM.py'

import optparse
import sys
import os
import numpy
from pylal import git_version
from pylal.fu_utils import getFOMdata,interferometers
disableGraphics=False
#Try importing GTK first if it fails import pylab non-interactive
#Try getting display env variable upon importing this module.
if os.getenv("DISPLAY") == None:
    #Non-interactive
    try:
        import matplotlib
        matplotlib.use("Agg")
        import pylab
    except Exception, errorInfo: #RuntimeError,ImportError:
        disableGraphics=True
        sys.stderr.write("Error trying to import NON-INTERACTIVE pylab!\n")
        sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
        sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
        sys.stderr.write("Pylab functionality unavailable!\n")
else:
    #Interactive
    try:
        import pylab
    except Exception, errorInfo: #RuntimeError,ImportError:
        disableGraphics=True
        sys.stderr.write("Error trying to import INTERACTIVE pylab!\n")
        sys.stderr.write("Exception Instance :%s\n"%(str(type(errorInfo))))
        sys.stderr.write("Exception Args     :%s\n"%(str(errorInfo.args)))
        sys.stderr.write("Pylab functionality unavailable!\n")

sys.path.append('@PYTHONLIBDIR@')
usage="""usage: %prog [options]"""

parser = optparse.OptionParser(usage,version=git_version.verbose_msg)
#Parse the command line to know what plots the user wants to make.
parser.add_option("-t","--gps-time",action="store",type="string",\
                  metavar="555444666",default=None,\
                  help="Specify the gps in INTEGER seconds. \
The resulting graphs will be plotted in such a way that t=0 \
will be the gps time given to this argument.")
parser.add_option("-w","--plot-windows",action="store",type="string",\
                  metavar="21600,7200",default="21600,7200",\
                  help="Specify this in seconds, using only positive \
numbers.  The default should be find for most people.")
parser.add_option("-i","--ifo-list",action="store",type="string",\
                  metavar="L1,H1",default=None,\
                  help="By default this plotting script will \
attempt to generate the FOMs for all plots and ifos it knows \
about see fu_utils.getFOMdata.channelDict for more details.")
parser.add_option("-g","--graphs",action="store",type="string",\
                  metavar="range,...,seis ",default=None,\
                  help="It is best to leave this one alone. If you \
want to explicitly specify the graphs to generate see \
fu_utils.getFOMdata.channelDict to figure out which graphs the code \
knows how to create.")
parser.add_option("-G","--graph-keys",action="store_true",default=False,\
                  help="Invoke this plot to display a list of graphs \
that can be generated.")
parser.add_option("-v","--verbose",action="store_true",
                  help="Set this flag to get some feedback during \
graph generation.")
parser.add_option("-p","--output-path",action="store",type="string",\
                  metavar="OUTPUT/PATH/FOR/FILES",default="./",\
                  help="Invoke this to specify where files are put.")
(opts,args) = parser.parse_args()
#
# Determine if given path is legit
#
outputPath=os.path.abspath(str(opts.output_path).strip())
if opts.verbose:
    sys.stdout.write("Figures to be saved to :%s"%outputPath)
    sys.stdout.flush()
if not os.path.exists(outputPath):
    if opts.verbose:
        sys.stdout.write("Creating path since it doesn't exist.\n")
        sys.stdout.flush()
    os.makedirs(outputPath)
    if not os.path.exists(outputPath):
        sys.stderr.write("Could not create needed path to save figures to.\n")
        sys.stderr.flush()
        os.abort()
#
# Setup the class to retrieve the data.
#
myFOM=getFOMdata(verbose=opts.verbose)
if opts.graph_keys:
    sys.stdout.write("Available graph keys :%s\n"%myFOM.getGraphKeys())
    sys.exit(0)
myFOM.setGPS(str(int(float(opts.gps_time))))
(preWindow,postWindow)=opts.plot_windows.strip().split(",")
myFOM.setWindows(preWindow,postWindow)
if opts.ifo_list==None:
    myIfos=interferometers
else:
    myIfos=opts.ifo_list.strip().split(",")
myFOM.setIfos(myIfos)
if opts.graphs==None:
    myGraphs=myFOM.getGraphList()
else:
    myGraphs=opts.graphs.strip().split(",")
myFOM.setGraphs(myGraphs)
#
# Retrieve the data
#
if opts.verbose:
    sys.stdout.write("Starting data fetching!\n")
graphData=myFOM.getData()
#
if opts.verbose:
    sys.stdout.write("Data retrieved.\n")
    for (key,obj) in graphData.iteritems():
        for (chan,data) in obj.iteritems():
            if data != None:
                sys.stdout.write("Graph: %s Trend: %s Length: %i\n"%(key,chan,len(data)))
                sys.stdout.flush()
            else:
                sys.stderr.write("Graph: %s Trend: %s Length: DataNotFound\n"%(key,chan))
                sys.stderr.flush()
#
# Generate the plots
#
#Meant as FOM_GPS_GRAPHNAME.png
filemask="FOM_%s_%s.png"
xRes=1600.0
xRes=1024
yRes=1200.0
yRes=768
myDPI=80.0
units={"h":float(3600.0),
       "s":float(1.0),
       "m":float(60.0),
       "d":float(86400.0)}
for thisGraph in myGraphs:
    if opts.verbose:
        sys.stdout.write("Creating graph : %s\n"%thisGraph)
    dataStreams=graphData[thisGraph].items()
    dataStreams.sort()
    myFigure=pylab.figure(figsize=(xRes/myDPI,yRes/myDPI),dpi=myDPI)
    pylab.title("%s t0 %s"%(thisGraph,opts.gps_time))
    myUnits="h"
    pylab.xlabel("Time (%s)"%myUnits)
    pylab.ylabel("Unknown")
    pylab.hold(True)
    for myLabel,data in dataStreams:
        if data != None:
            metaDataDict=data.metadata.todict()
            dt=float(metaDataDict["dt"])
            dataSegments=metaDataDict["segments"]
            myData=data.tolist()
            myStart=float(dataSegments[0][0])
            myStop=myStart+len(myData)*float(dt)
            timeVector=numpy.arange(myStart,myStop,dt)
            #Convert the time units
            timeVector=[(float(x)-float(opts.gps_time))/units[myUnits] \
                        for x in timeVector]
            if len(timeVector) != len(myData):
                sys.stderr.write(\
                    "Inconsistent Vector lengths Channel %s: %s vs %s\n"\
                    %(myLabel,len(timeVector),len(myData)))
            pylab.plot(timeVector,myData,label="%s"%myLabel)
    imageFilename=os.path.normpath(outputPath+"/"+filemask%(opts.gps_time,thisGraph.replace(" ","-")))
    #Create thumbnail without legend first!
    myScale=0.33
    myFigure.set_size_inches(myScale*(xRes/myDPI),myScale*(yRes/myDPI))
    myFigure.set_dpi(myScale*myDPI)
    pylab.savefig(imageFilename.replace(".png",".thumb.png"))
    #If inspiral plot put legend lower left
    if thisGraph.lower().__contains__("inspiral"):
        pylab.legend(loc=3)
    else:
        pylab.legend()
    pylab.grid(True)
    #Save large image
    myFigure.set_size_inches(1.0*(xRes/myDPI),1.0*(yRes/myDPI))
    myFigure.set_dpi(1.0*myDPI)
    pylab.savefig(imageFilename)    
    pylab.close()
#
# Done
#
