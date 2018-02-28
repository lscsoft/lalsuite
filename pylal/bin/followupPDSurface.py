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
__author__ = "Cristina Valeria Torres <cristina.torres@ligo.org>"


"""
This script lets you specify a channel and gps markers, from this a
ascii file is generated and written to disk.  This is a simple tool
to allow novice persons to analyse data without understanding the
complete LIGO data handling infrastructure.
"""
from StringIO import StringIO
from commands import getstatusoutput
from glue import lal
from numpy import mean,std,median,min,max,\
     abs,array,log,floor,power,arange,\
     floor,ceil,searchsorted,empty,isnan,\
     interp,zeros,arctan,pi,ones,size
from numpy.ma import masked_array
from pylal import frutils
from sys import stderr,stdout,exit
import optparse
import matplotlib
matplotlib.use("Agg")
import pylab
from scipy import interp
from pylal import followup_utils

#
# Misc interval methods
#
def cutStreamData(tData,yData,t0,t1):
    """
    Returns two vectors (tVal,yVal) as
    a tuple. Returns (None,None) if cuts faile.
    """
    if len(tData) != len(yData):
        print "Warning data vectors do not match!"
        print "tData %i: yData %i"%(len(tData),len(yData))
        print "Rescaling to size of time vector assumed correct!."
        yData=yData[0:len(tData)]
    cut0=searchsorted(tData,t0,side='left')
    cut1=searchsorted(tData,t1,side='right')
    return (tData[cut0:cut1],yData[cut0:cut1])

def getFrObjectAsLists(FrDataStream):
    """
    Returns two lists as (tStamps,dataPoints)
    """
    metaDataDict=FrDataStream.metadata.todict()
    dT=float(metaDataDict["dt"])
    segmentInfo=metaDataDict["segments"]
    dataPoints=FrDataStream.tolist()
    myStart=segmentInfo[0][0]
    myStop=myStart+len(dataPoints)*dT
    tStamps=arange(myStart,myStop,dT)
    return (tStamps,dataPoints)
#End DEF

def writeStream(FrDataStream,filename):
    """
    Writes out the dBata stream to a text file
    """
    fp=open(filename,'w')
    tmpMetaData=FrDataStream.metadata.todict()
    for myKey,myVal in tmpMetaData.iteritems():
        fp.write("#%s:%s\n"%(myKey,myVal))
    #
    # Write out data
    #
    for datum in FrDataStream.tolist():
            fp.write("%s\n"%datum)
    fp.close()

def cutStream(FrDataStream,t0,t1):
    """
    Uses the metadata and cuts the frame stream.  It returns
    a numpy array to you.  The meta data is adjusted accordingly
    in the structure.
    """
    tData=getTimeStamps(FrDataStream)
    cut0=searchsorted(tData,t0,side='left')
    cut1=searchsorted(tData,t1,side='right')
    cutStream=FrDataStream[cut0:cut1]
    return cutStream

def getTimeStamps(FrDataStream):
    """
    Returns a list of time stamps tabulated from
    meta data in input.
    """
    metaDataDict=FrDataStream.metadata.todict()
    dT=float(metaDataDict["dt"])
    segmentInfo=metaDataDict["segments"]
    myStart=segmentInfo[0][0]
    myStop=myStart+len(FrDataStream)*dT
    return arange(myStart,myStop,dT)

def build3DHistogramArray(inputA=None,xBinVector=None,yBinVector=None):
    """
    This builds and bins up the DCPD exposure surface.
    The input is a time parameterized array
    [beamlist,beampitchlist,beamyawlist]
    """
    threeDData=list()
    for i,aCol0 in enumerate(inputA[0]):
        threeDData.append((inputA[0][i],\
                           inputA[1][i],\
                           inputA[2][i]))
    #Sort the input structure by col0
    threeDData.sort()
    zMapLists=empty([len(xBinVector),len(yBinVector)],dtype=object)
    for ii in range(0,len(xBinVector)):
        for jj in range(0,len(yBinVector)):
            zMapLists[ii][jj]=list()
    for xIndex in range(0,len(xBinVector)-1):
        xLow=xBinVector[xIndex]
        xHigh=xBinVector[xIndex+1]
        xData=[a for a,b,c in threeDData]
        # Find all points that fit this X column
        dataSubset=threeDData[searchsorted(xData,xLow,side='left'):\
                              searchsorted(xData,xHigh,side='right')]
        yDataSubset=[(b,c) for a,b,c in dataSubset]
        yDataSubset.sort()
        yData=[b for b,c in yDataSubset]
        for yIndex in range(0,len(yBinVector)-1):
            yLow=yBinVector[yIndex]
            yHigh=yBinVector[yIndex+1]
            dataBinMatch=yDataSubset[searchsorted(yData,yLow,side='left'):\
                                     searchsorted(yData,yHigh,side='right')]
            zDataMatch=[c for b,c in dataBinMatch]
            zMapLists[xIndex][yIndex].extend(zDataMatch)
    return zMapLists

def scatterPointer(A=None,B=None):
    """
    Expects two ordered pairs as tuples (x1,y1) and (x2,y2)...
    Then this defines an angle.
    From this angle we orient a triangle point and return this as
    a tuple of requested size.  This can be plotted by the scatter
    plot to point in a particular orientation.
    The direction right(0,0)->(1,0) is angle 0 radians, then we rotate
    counter clockwise from there... What is returned in a three
    element tuple (numsides,style,angle) which can be put into
    a plot call via marker=X as a **kwarg
    """
    if A == None or B == None:
        return (3,0,0)
    if( A != type(tuple()) and len(A) != 2) \
        or \
      ( B != type(tuple()) and len(B) != 2):
        return (3,0,0)
    #
    # Calculate orientation of triangle
    # Ang = arcsin(dY/dX)
    dY=float(B[-1]-A[-1])
    dX=float(B[0]-A[0])
    if dY>=0 and dX>0:
        myAngle=arctan(dY/dX)
    elif dY>0 and dX<0:
        myAngle=pi-arctan(dY/dX)
    elif dY<0 and dX<0:
        myAngle=arctan(dY/dX)+pi
    elif dY<0 and dX>0:
        myAngle=(2.0*pi)-arctan(dY/dX)
    elif dX==0 and dY>=0:
        myAngle=(pi/2.0)
    elif dX==0 and dY<0:
        myAngle=(3.0/2.0)*pi
    elif dX>0 and dY==0:
        myAngle=0
    elif dX<0 and dY==0:
        myAngle=pi
    else:
        myAngle=0
    return (3,0,myAngle)
        
def convert2DHistTo3ColVectors(myMatrix,xBins=None,yBins=None):
    """
    Create three vectors X,Y,Z for ease of use in functs like contour,
    and scatter instead of imshow!
    If [x,y]Bins==None use bin index for X, Y
    Input matrix assumed square!!
    """
    if xBins==None:
        xBins=range(len(myMatrix))
    if yBins==None:
        yBins=range(len(myMatrix[0]))
    if len(xBins)*len(yBins) != myMatrix.size:
        raise Exception, "Input matrix and bin vectors disagree!"
    tX=empty(myMatrix.size,type(xBins[0]))
    tY=empty(myMatrix.size,type(yBins[0]))
    tZ=empty(myMatrix.size,type(myMatrix[0][0]))
    myIndex=0
    for ii,iVal in enumerate(xBins):
        for jj,jVal in enumerate(yBins):
            tX[myIndex]=iVal
            tY[myIndex]=jVal
            tZ[myIndex]=myMatrix[ii][jj]
            myIndex=myIndex+1
    return tX,tY,tZ

usage = """usage: %prog [options]"""
validIFOList=["L","H","V"]
parser = optparse.OptionParser(usage)
#Specify all needed options
parser.add_option("-l","--laser-channel",action="store",\
                  type="string",default="dummy",\
                  help="Name of the laser channel.")
parser.add_option("-p","--laser-pitch",action="store",\
                  type="string",default=None,\
                  help="Name of the laser pitch channel.")
parser.add_option("-y","--laser-yaw",action="store",\
                  type="string",default=None,\
                  help="Name of the laser yaw channel.")
parser.add_option("-s","--gps-t0",action="store",type="string",\
                   metavar="gpsT0",default=None,\
                   help="Specify the GPS t0 time of the \
for data interval of interest")
parser.add_option("-e","--gps-window",action="store",type="string",\
                   metavar="gpsWindow",default="0.5,0.5",\
                   help="Specify the integer window time pair \
(InFront,Behind) for the data point t0 of interest. Default is \
pm 0.5,0.5 second(s)")
parser.add_option("-d","--gps-history",action="store",type="string",\
                   metavar="gpsHistory",default="300,300",\
                   help="Specify a second window pair \
(InFront,Behind)in integer seconds for constructing the \
PSD 2D histograms to determine surface senstivity for DCPDs. \
Default is 300,300 seconds")
parser.add_option("-i","--ifo",action="store",type="string",\
                   metavar="IFO",default=None,help="Specify the \
IFO of interest as one of the following from the set of valid IFO \
args %s"%validIFOList)
parser.add_option("-t","--type",action="store",type="string",\
                   metavar="FRAMETYPE",default="R",\
                   help="Specify the type of frame to make a query \
for.  See the ligo_data_find help for a complete list of \
frame types. Default to 'R' raw frames.")
parser.add_option("-u","--url-type",action="store",type="string",\
                   metavar="URLTYPE",default="file",\
                   help="The default URL type to query for is local \
file.  If you wish to change it to another availble type see \
the ligo_data_find help for the available types.")
parser.add_option("-o","--output-filename-base",action="store",\
                   type="string",metavar="NAME",default="QuickCheck",\
                   help="Set this flag to override the automagic \
naming of files to something other than channel_gpsstart_deltaT.ascii\
for the raw data and channel_gpsT0_deltaT.png for the plots.\
The specified name here will always get either txt or png as extension.")
parser.add_option("-b","--bin-count",action="store",type="string",\
                  metavar="150",default="150",\
                  help="Specify bin count of smallest axis.")
parser.add_option("-D","--dumps-off",action="store_true",\
                  default=False,\
                  help="Set this flag to disable creation of \
raw data ascii files.  This will speed things up.")
parser.add_option("-v","--verbose",action="store_true",\
                   default=False,\
                   help="Set this flag to make the software \
output progress information to the shell.")
#
# Process command arguments
#
(opts,args) = parser.parse_args()
if (opts.url_type == None) or \
   (opts.type == None) or \
   (opts.gps_t0 == None) or \
   (opts.laser_channel==None) or \
   (opts.laser_pitch==None) or \
   (opts.laser_yaw==None):
    stderr.write("Invalid arguments passed or insufficient \
    arguments given!\n")
    exit(1)
#
# Create data spigot
#
(gps_frontHistory,gps_backHistory)=opts.gps_history.split(",")
gpsStart=floor(float(opts.gps_t0)-float(gps_frontHistory))
gpsEnd=ceil(float(opts.gps_t0)+float(gps_backHistory))
(gps_frontWindow,gps_backWindow)=opts.gps_window.split(",")
gpsA=floor(float(opts.gps_t0)-float(gps_frontWindow))
gpsB=ceil(float(opts.gps_t0)+float(gps_backWindow))
gpsT0=float(opts.gps_t0)
beamName=str(opts.laser_channel)
beamPitch=str(opts.laser_pitch)
beamYaw=str(opts.laser_yaw)
if opts.verbose:
    stdout.write("Accessing and Reading in datasets.\n")
#
# Create connection to frame cache of data
#
beamSpigot=followup_utils.connectToFrameData(gpsStart,
                                             gpsEnd,
                                             opts.type,
                                             opts.ifo,
                                             opts.url_type,
                                             True,
                                             opts.verbose)
#
# Grab all data for beam intensity, beam pitch and beam yaw 
#
origData=dict()    
origData["pitch"]=beamSpigot.getDataStream(beamPitch,gpsStart,gpsEnd)
origData["yaw"]=beamSpigot.getDataStream(beamYaw,gpsStart,gpsEnd)
if beamName=="dummy":
    # In the case of needing dummy data, we will just load up the pitch
    # data vector and replace the numpy data structure with ones!
    origData["beam"]=frutils.TimeSeries(ones(size(origData["pitch"])),origData["pitch"].metadata)
else:
    origData["beam"]=beamSpigot.getDataStream(beamName,gpsStart,gpsEnd)
#
# Tabulate the bins for binning the PD surface
#
binCount=int(opts.bin_count)
minRes=min([origData["pitch"].max()-origData["pitch"].min(),\
            origData["yaw"].max()-origData["yaw"].min()])/binCount
## minValue=min([origData["yaw"].min(),origData["pitch"].max()])
## maxValue=min([origData["yaw"].max(),origData["pitch"].max()])
## bBox=(minValue,maxValue,minValue,maxValue)
#[Lt,Rt,Bot,Top]
bBox=(origData["yaw"].min(),origData["yaw"].max(),\
      origData["pitch"].min(),origData["pitch"].max())
yawBinVector=pylab.arange(bBox[0],bBox[1],minRes)    
pitchBinVector=pylab.arange(bBox[2],bBox[3],minRes)
#
# Retabulate the data streams with(lowest df) a common time sampling max(dT)
#
minDT=min([x["dt"] for x in \
           [origData["beam"].metadata.todict(),\
            origData["pitch"].metadata.todict(),\
            origData["yaw"].metadata.todict()]])
tPrime=arange(gpsStart,gpsEnd,minDT)
#
# Cast all revised data into a dict object
#
myLabel=dict()
myLabel["time"]="Time(s)"
myLabel["beam"]=beamName
myLabel["pitch"]=beamPitch
myLabel["yaw"]=beamYaw
myData=dict()
myData["time"]=tPrime
for myKey in origData.keys():
    if myKey not in "time":
        myData[myKey]=interp(myData["time"],\
                             getTimeStamps(origData[myKey]),\
                             origData[myKey])
if opts.verbose:
    stdout.write("Creating associated data to plot.\n")
#
# Create 2D histogram of beam intensities
#
beamIntensityHistData=build3DHistogramArray([myData["yaw"],
                                             myData["pitch"],
                                             myData["beam"]],
                                            yawBinVector,
                                            pitchBinVector)
#
# Convert to histogram of beam location
#
rangeLocationData=list()
beamLocationHistData=empty(beamIntensityHistData.shape,\
                           dtype=float)
for ii in range(beamIntensityHistData.shape[0]):
    for jj in range(beamIntensityHistData.shape[1]):
        beamLocationHistData[ii][jj]=minDT*float(len(beamIntensityHistData[ii][jj]))
        if not isnan(beamLocationHistData[ii][jj]):
            rangeLocationData.append(beamLocationHistData[ii][jj])
#
# Convert to histogram of median beam power per coord
#
beamMedianHistData=empty(beamIntensityHistData.shape,\
                         dtype=float)
rangeMedianData=list()
for ii in range(beamIntensityHistData.shape[0]):
    for jj in range(beamIntensityHistData.shape[1]):
        beamMedianHistData[ii][jj]=median(beamIntensityHistData[ii][jj])
        if not isnan(beamMedianHistData[ii][jj]):
                     rangeMedianData.append(beamMedianHistData[ii][jj])
#
# Create a mapping of beam motion for smaller time interval
#
mySnipData=dict()
mySnipData["time"]=arange(gpsA,gpsB,minDT)
mySnipData["time"]=arange(float(opts.gps_t0)-float(gps_frontWindow),\
                          float(opts.gps_t0)+float(gps_backWindow),
                          minDT)
for myKey in myData.keys():
    if "time" not in myKey:
        if beamName == "dummy" and myKey == "beam":
            mySnipData[myKey]=ones(size(mySnipData["time"]))
        else:
            tmpData=beamSpigot.getDataStream(myLabel[myKey],gpsA,gpsB)
            mySnipData[myKey]=interp(mySnipData["time"],
                                     getTimeStamps(tmpData),
                                     tmpData)
#
# Setup figure
#
xRes=1600.0
yRes=1200.0
myDPI=80.0
myHandles=list()
myFigure=pylab.figure(figsize=(xRes/myDPI,yRes/myDPI),dpi=myDPI)
#
# Figure of PD beam time on
#
mySubPlot=pylab.subplot(1,3,1)
# Setup Mask for "not paint" zero time bins
myColorMap=matplotlib.cm.jet
myColorMap.set_bad('w',0.0)
maskedBeamLocationHistData=masked_array(beamLocationHistData,
                                        beamLocationHistData<=0.0)
myHandles.append(pylab.imshow(maskedBeamLocationHistData.transpose(),\
                              interpolation='nearest',\
                              extent=bBox,\
                              aspect='auto',\
                              cmap=myColorMap,\
                              origin='lower'))
myYaw=interp(gpsT0,myData["time"],myData["yaw"])
myPitch=interp(gpsT0,myData["time"],myData["pitch"])
myIntensity=interp(gpsT0,myData["time"],myData["beam"])
myYawPlusDT=interp(gpsT0+minDT,myData["time"],myData["yaw"])
myPitchPlusDT=interp(gpsT0+minDT,myData["time"],myData["pitch"])

starSize=75
myBeam=interp(gpsT0,myData["time"],myData["beam"])
mySubPlot.annotate('*',xy=(myYaw,myPitch),xycoords='data',\
                   size=starSize,color='red')
pylab.title("Beam Time on PD (%s,%s,%s)"%(gpsT0,float(gps_frontHistory),float(gps_backHistory)))
#Need to determine colorbar range
CB1=pylab.colorbar(orientation='horizontal')
CB1.set_label("Time(s)")
pylab.clim(min(rangeLocationData),max(rangeLocationData))
pylab.xlabel(myLabel["yaw"])
pylab.ylabel(myLabel["pitch"])
#
# Figure of PD intensity mapping
#
mySubPlot2=pylab.subplot(1,3,2)
myHandles.append(pylab.imshow(beamMedianHistData.transpose(),\
                              interpolation='nearest',\
                              extent=bBox,\
                              aspect='auto',\
                              origin='lower'))
mySubPlot2.annotate('*',xy=(myYaw,myPitch),xycoords='data',\
                   size=starSize,color='red')
prevFig=myHandles[-1]
pylab.title("Median Beam Intensity on PD (%s,%s,%s)"%(gpsT0,float(gps_frontHistory),float(gps_backHistory)))
#Need to determine colorbar range
CB2=pylab.colorbar(orientation='horizontal')
CB2.set_label("Counts")
pylab.clim(min(rangeMedianData),max(rangeMedianData))
pylab.xlabel(myLabel["yaw"])
pylab.ylabel(myLabel["pitch"])
#
# Figure a beam trace of short interval around t0
#
#
# Figure of PD intensity mapping
#
mySubPlot3=pylab.subplot(1,3,3)
myHandles.append(pylab.scatter(mySnipData["yaw"],
                               mySnipData["pitch"],
                               c=mySnipData["beam"],
                               edgecolors=None,
                               linewidth=0,
                               hold=True,\
                               label="Beam Path"))
CB3=pylab.colorbar(orientation='horizontal')
pylab.clim(min(rangeMedianData),max(rangeMedianData))
CB3.set_label("Counts")
symbolSize=250
myMarkerT0=scatterPointer((myYaw,myYawPlusDT),(myPitch,myPitchPlusDT))
print "Coord :",myYaw,myYawPlusDT,myPitch,myPitchPlusDT,myMarkerT0
myHandles.append(pylab.scatter([myYaw],\
                               [myPitch],\
                               s=symbolSize,\
                               c='m',\
                               marker=myMarkerT0,
                               linewidth=1,\
                               label="t0"))
myMarkerStart=scatterPointer((mySnipData["yaw"][0],mySnipData["pitch"][0]),\
                             (mySnipData["yaw"][1],mySnipData["pitch"][1]))                        
myHandles.append(pylab.scatter([mySnipData["yaw"][0]],\
                               [mySnipData["pitch"][0]],\
                               s=symbolSize,\
                               c='g',\
                               marker=myMarkerStart,\
                               linewidth=1,\
                               label="Start"))
myMarkerStop=scatterPointer((mySnipData["yaw"][-1],mySnipData["pitch"][-1]),\
                        (mySnipData["yaw"][-2],mySnipData["pitch"][-2]))
myHandles.append(pylab.scatter([mySnipData["yaw"][-1]],\
                               [mySnipData["pitch"][-1]],\
                               s=symbolSize,\
                               c='r',\
                               marker=myMarkerStop,\
                               linewidth=1,\
                               label="Stop"))                               
pylab.title("Beam trace on PD face at %s,%s,%s"%(gpsT0,float(gps_frontWindow),float(gps_backWindow)))
#Need to determine colorbar range
pylab.xlabel(myLabel["yaw"])
pylab.ylabel(myLabel["pitch"])
pylab.legend()
#
# Print out beam location to screen
#
print "Beam location at %s is (%s,%s) with intensity measure %s"%\
      (gpsT0,\
       myYaw,\
       myPitch,\
       myIntensity)
print "Beam location at %s is (%s,%s) with intensity measure %s"%\
      (mySnipData["time"][0]-gpsT0,\
       mySnipData["yaw"][0],\
       mySnipData["pitch"][0],\
       mySnipData["beam"][0])
print "Beam location at %s is (%s,%s) with intensity measure %s"%\
      (gpsT0-mySnipData["time"][-1],\
       mySnipData["yaw"][-1],\
       mySnipData["pitch"][-1],\
       mySnipData["beam"][-1])

#
# Save the raw data and plots to disk.
#
if not opts.dumps_off:
    if opts.verbose:
        stdout.write("Writing original(unresampled) data as ascii files to disk.\n")
    for myKey in origData.keys():
        if "time" not in myKey:
            myFilename="%s_%s_%s_%s.ascii"%(opts.output_filename_base,
                                            myLabel[myKey],
                                            gpsT0,
                                            opts.gps_history)
        writeStream(origData[myKey],myFilename)
#
# Save the figure to the disk
#
if opts.verbose:
    stdout.write("Saving graphics to disk.\n")
myFilename="%s_%s_%s_%s.png"%(opts.output_filename_base,
                              myLabel["beam"],
                              gpsT0,
                              opts.gps_history)
pylab.savefig(myFilename)
pylab.close()

if opts.verbose:
    stdout.write("Graphics prepared!\n")

exit(0)
