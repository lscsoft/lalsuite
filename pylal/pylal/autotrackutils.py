"""
  * Copyright (C) 2004, 2005 Cristina V. Torres
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with with program; see the file COPYING. If not, write to the
  *  Free Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,
  *  MA  02111-1307  USA
"""
__author__ = 'Cristina Torres <cristina.torres@ligo.org>'

import os
import numpy
import copy
import sys
try:
  import sqlite3 as sqlite
except ImportError:
  import sqlite
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

"""
This file is intended to provide the autotrack utilities required to
process the results of the tracksearch hybrid MDGC codes.
"""

def writeTGNListToFile(filename="default",tgnList=None):
  """
  This writes the information in the tgnList to a simple text
  file. Actually we write two files.  filename.TGN and
  filename.AUX.  In the AUX file is information like time stamps,
  and other misc information from the TGNs that are dumped to disk!
  """
  if tgnList==None:
    return
  #Append extensions AUX and TGN
  fp0=open(filename+".TGN",'wt')
  fp1=open(filename+".AUX",'wt')
  for tgn in tgnList:
    a,b=tgn.reverseCreateFromString()
    id=tgn.getID()
    count=tgn.getMemberCount()
    vol=tgn.getVolume()
    fp0.write("%s\n%s\n"%(a,b))
    fp1.write("%s\t%s\t%s\n"%(id,count,vol))
  fp0.close()
  fp1.close()
  #End writeTGNListToFile()
  
def generateTGNListFromFile(filename="",timestamp=str("0")):
  """
  This method opens a text file containing multiple TGNs and creates
  a list of TGN objects.  It is assumed that this list of TGNs was
  created at the same time for each TGN in the file.
  """
  fp=open(filename,'rt')
  rawData=fp.readlines()
  fp.close()
  linesPerTGN=2
  if rawData.__len__().__mod__(linesPerTGN) != 0:
      raise autotrackError("File appears inconsistent!")
  eCount=rawData.__len__().__div__(linesPerTGN)
  TGNList=list()
  for index in range(0,eCount):
      thisTGN=TGN()
      thisTGN.createFromString(rawData[index*2],rawData[index*2+1])
      thisTGN.__setBirthDate__(timestamp)
      thisTGN.__setGPS__(int(timestamp))
      TGNList.append(thisTGN)
  return TGNList
    #End

def getQualityTable(backEndName="tracksearch"):
  """
  Method is a complex case statement that selects out given a
  back-end name the properties that the TGN should be capable of
  tracking. Depends on the output formatting of tracksearchMDGCv2.m
  file
  """
  if backEndName.strip().lower() == "tracksearch":
      #Fetch the property fields method out of the tracksearchutils module
      from lalapps.tracksearchutils import candidateListProperties
      tmpProperties=candidateListProperties()
      #We assume that the below elements to be defined were
      #processed in tracksearchMDGCv2 et all files and are in that
      #order left to right columns. The column values specified here
      #are those needed to access any particular data field from the
      #
      #Structure ["tag",TGNarrayIndex]
      propertiesToKeep=[
          ["pp",0],
          ["y",1],
          ["b",2],
          ["d",3],
          ["a",4],
          ["r",5],
          ["rr",6],
          ["w",7],
          ["e",8],
          ["ww",9],
          ["ee",10]
          ]
      propertiesHandDefined=[
          ["ht","Time symmetry",11],
          ["hf","Freq symmetry",12]
          ]
      #Create a 3 list with text labels and array locations!
      finalPropertyList=[]
      for myRow in propertiesToKeep:
          for thatRow in tmpProperties:
              if thatRow[0].lower().strip() == \
                      myRow[0].lower().strip():
                  finalPropertyList.append([myRow[0],thatRow[1],myRow[1]])
      for myRow in propertiesHandDefined:
          finalPropertyList.append(myRow)
  return finalPropertyList
#End getQualityTable()

def coalescetuplelist(inputList=list()):
    """
    Takes an input list of tuples representing start,stop and merges
    them placing them in time order when returning the coalesced list of
    tuples.
    """
    outputList=list()
    if type(inputList) != type(list()):
      sys.stderr.write("Wrong variable type passed as argument in\
 autotrackutils.__coalescetuplelist__()\n")
      return None
    if inputList.__len__() < 1:
      return  inputList
    inputList.sort()
    while inputList:
        segA=inputList.pop()
        overlap=True
        #Assume next segment overlaps segA
        while overlap:
            #Pop of next segment if available
            if inputList.__len__() > 0:
                segB=inputList.pop()
            else:
                #No overlap possible no segs left!
                segB=(-1,-1)
                overlap=False
            #Three cases of intersection
            #Overlap Left
            if (
                (segB[0]<= segA[0] <= segB[1])
                and
                (segA[1] >= segB[1])
                ):
                segA=(segB[0],segA[1])
            #Overlap Right
            elif (
                  (segB[0]<= segA[1] <= segB[1])
                  and
                  (segA[1] <= segB[0])
                 ):
                segA=(segA[0],segB[1])
            #Bridge over
            elif (
                (segB[0]<=segA[0])
                and
                (segB[1]>=segA[1])
                ):
                segA=(segB[0],segB[1])
            else:
                #Put segment back there was no overlap!
                if not((-1,-1)==segB):
                  inputList.append(segB)
                overlap=False
        outputList.append(segA)
        outputList.sort()
    return outputList
#End def coalescetuplelist(input=list()):

def activeThreadIntervals(tgnList=list(),tgnSize=int(300)):
  """
  A simple coalesce feature useful to autotrack to return a list of
  segments which define intervals that the thread was active.
  """
  tupleList=list()
  for tgn in tgnList:
    start=float(tgn.getBirthDateText())
    stop=start+tgnSize
    tupleList.append((start,stop))
  activeListofTuples=coalescetuplelist(tupleList)
  return activeListofTuples
#End activeThreadIntervals(tgnList=list(),tgnSize=int(300)):

class autotrackError(Exception):
  """
  Generic class raise this flag if we do not have a better
  flag name.
  """
  #End class

class autotrackDefineMismatchError(Exception):
  """
  Custom error exceptions for these class objects.
  """
  #End class autotrackDefineError

class autotrackDefineIdMismatchError(Exception):
  """
  Custom error handling of mismatched ID
  """
  #End class autotrackDefineIdMismatchError


class plotTGN:
  """
  This class is a wrapper for plotting TGNs.  It can plot TGNs on a
  2D plot of requested parameter options.
  """
  def __init__(self,iTGN=None,thread=bool(False)):
    """
    Ths method should be invoked blind with either a single TGN or
    a list of TGNS.  Invoking without a TGN is still possible.  A
    second argument lets the plotting code know if we are
    comparing different TGNs or plotting a thread of TGNs.
    """
    self.isThread=thread
    if (iTGN != None and type(iTGN) != type([])):
        raise autotrackError("plotTGN expects a list of TGN objects.")
    self.tgnList=iTGN
    #Default quantities to plot
    self.defaultX="d"
    self.defaultY="y"
    #Default member scale for plots 1px == 100members
    self.memberSteps=100;
  #End __init__()

  def setXQuantity(self,newX="d"):
    """
    Using the single character string see the output of
    getQualityTable() select the quantity to show on the Xaxis.
    """
    self.defaultX=newX
  #End setXQuantity

  def setYQuantity(self,newY="y"):
    """
    Using the single character string see the output of
    getQualityTable() select the quantity to show on the Yaxis.
    """
    self.defaultY=newY
  #End setYQuantity

  def createFigure(self,filename=None):
    """
    This function provides an interface to python plotting that
    will plot a TGN or thread of TGNs given to the plotTGN class.
    """
    if self.tgnList == None:
        raise autotrackError("Figure creation requested, no TGNs defined.")
    pylab.figure()
    if not self.isThread:
        self.__primativePlotTGNs__()
    if self.isThread:
        self.__primativePlotThread__()
    if filename == None:
      pylab.show()
    else:
      pylab.savefig(filename)
      pylab.close()
  #End createFigure()

  def createReportFigures(self,filename=None):
    """
    This function creates a single plot with subplot of all know
    properties plotted against each other to see if there are 
    any relationships between TGNs
    """
    if self.tgnList == None:
      raise autotrackError("Figure creation requested, no TGNs defined.")
    figureCount=0
    if filename != None:
      [file,ext]=os.path.splitext(filename)
    #Get property list
    shortName=[element[0] for element in getQualityTable()]
    plotCount=shortName.__len__()*shortName.__len__()
    #Row image count
    countRowImages=4
    #Col image count
    countColImages=int(shortName.__len__()).__div__(countRowImages)
    if int(shortName.__len__()).__mod__(countRowImages) > 0:
      countColImages=countColImages+1
    for myXAxis in shortName:
      myHandle=pylab.figure()
      countImage=1
      for myYAxis in shortName:
        if myXAxis != myYAxis:
          self.setXQuantity(myXAxis)
          self.setYQuantity(myYAxis)
          pylab.subplot(countRowImages,countColImages,countImage)
          countImage=countImage+1
          self.__primativePlotTGNs__(True)
      if filename != None:
        newFilename="%s_%s%s"%(file,str(figureCount).zfill(3),ext)
        pylab.savefig(newFilename,dpi=160)
        pylab.close()
      figureCount=figureCount+1
    if filename == None:
      pylab.show()
  #End createReportFigures()

  def __getIndexAndLabels__(self):
    """
    Returns index and labels for X and Y axis
    """
    for elem in getQualityTable():
        if elem[0].strip().lower() == self.defaultX:
            xIndex=elem[2]
            xLabel=elem[1]
        if elem[0].strip().lower() == self.defaultY:
            yIndex=elem[2]
            yLabel=elem[1]
    return [xIndex,xLabel,yIndex,yLabel]
  #End __getIndexAndLabels__()

  def __primativePlotTGNs__(self,bare=bool(False)):
    """
    Is a macro of plotting commands that takes a list of TGNs that
    plots each of these individually as a collection of points.
    Creates a figure plotting the thread of list of TGNs using the
    centroid and an X,Y error bars.  Take a optional boolean to
    make the plot not include a title and legend
    """
    #Determine the index that corresponds to X and Y quantities
    xIndex=0
    yIndex=0
    xLabel="NULL"
    yLabel="NULL"
    (xIndex,xLabel,yIndex,yLabel)=self.__getIndexAndLabels__()
    plotValues=list()
    gpsTimesInList=list()
    for thisTGN in self.tgnList:
        label=str(thisTGN.getID())
        #Get the X,Y property
        (xC,xE)=thisTGN.getCentroidErrorViaIndex(xIndex)
        (yC,yE)=thisTGN.getCentroidErrorViaIndex(yIndex)
        plotValues.append([xC,yC,xE,yE,label])
        gpsTimesInList.append(thisTGN.getGPS())
    for x,y,ex,ey,txtLabel in plotValues:
        pylab.errorbar(x,y,xerr=ex,yerr=ey,label=txtLabel,marker='o')
    pylab.xlabel(str(xLabel))
    pylab.ylabel(str(yLabel))
    if not bare:
      pylab.title("TGNs: %i"%(min(gpsTimesInList)))
      pylab.legend()
  #End __primativePlotTGNs__

  def __primativePlotThread__(self,bare=bool(False)):
    """
    Is a macro of plotting commands that takes a list of TGNs that
    plots each of these individually as a collection of points.
    Creates a figure plotting the thread of list of TGNs using the
    centroid and an X,Y error bars
    """
    #Determine the index that corresponds to X and Y quantities
    xIndex=0
    yIndex=0
    xLabel="NULL"
    yLabel="NULL"
    for elem in getQualityTable():
        if elem[0].strip().lower() == self.defaultX:
            xIndex=elem[2]
            xLabel=elem[1]
        if elem[0].strip().lower() == self.defaultY:
            yIndex=elem[2]
            yLabel=elem[1]
    plotValues=list()
    gpsTimesInList=list()
    for thisTGN in self.tgnList:
        txtLabel=str(thisTGN.getGPS())
        #Get the X,Y property
        (xC,xE)=thisTGN.getCentroidErrorViaIndex(xIndex)
        (yC,yE)=thisTGN.getCentroidErrorViaIndex(yIndex)
        plotValues.append([xC,yC,xE,yE,txtLabel])
        gpsTimesInList.append(thisTGN.getGPS())
    xVec=list()
    yVec=list()
    for x,y,ex,ey,txtLabel in plotValues:
        xVec.append(x)
        yVec.append(y)
    pylab.plot(xVec,yVec,label="Time Line")
    for x,y,ex,ey,txtLabel in plotValues:
        pylab.errorbar(x,y,xerr=ex,yerr=ey,label=txtLabel,marker='o')
    pylab.xlabel(str(xLabel))
    pylab.ylabel(str(yLabel))
    if not bare:
      pylab.title("TGN Thread %i"%(min(gpsTimesInList)))
      pylab.legend()
  #End __primativePlotThread__


class TGN:
  """
  This class provides the definition of a single defined autotrack
  defined neighborhood.  We use these neighborhoods to track
  how the instrument behavior groups change, appear or disappear.
  """
  def __init__(self):
      """ This method initializes an empty neighborhood.
      To populate the neightborhood invoke the proper method
      depending on the source of the data, ie ascii mysql etc
      """
      self.delim=" "
      self.idNum=float(-1)
      self.density=float(0)
      self.memberCount=0
      self.volume=0
      self.colCount=0
      self.birthdate=str("-1")
      self.gpsStamp=int(-1)
      self.discoverer=str()
      self.lastSeen=float(0)
      self.center=None
      self.bound=None
      #Convention below should match tracksearchUtils conventions
      #Only keeping the lines that apply from method glitchDB
      #variable named glitchDatabaseEntry method should be expanded
      #to not be hardwired to backend defintions
      self.qualities=getQualityTable("tracksearch")
  #End __init__ method

  def __setID__(self,inputArg):
      """
      Should set the ID numeric field
      """
      if type(float(0)) != type(inputArg):
          inputArg=float(inputArg)
      self.idNum=inputArg
  #End

  def getID(self):
      """
      Fetches the numeric ID assigned to this TGN instance
      """
      return self.idNum
  #End

  def __setBirthDate__(self,bDate=str("")):
      """
      Set a text string which is the GPS birthdate,
      as closely as possible for this neighborhood.
      """
      if type(str()) != type(bDate):
          bDate=str(bDate)
      self.birthdate=bDate
  #End

  def __setGPS__(self,gps=int(0)):
      """
      Set a text string which is the GPS birthdate,
      as closely as possible for this neighborhood.
      """
      if type(int()) != type(gps):
          raise autotrackError("GPS time not a INT type.")
      self.gpsStamp=gps
  #End

  def getBirthDateText(self):
      """
      This retrieves the text string birthdate stored in TGN
      instance.
      """
      return self.birthdate
  #End

  def getGPS(self):
      """
      This retrieves the integer GPS birthdate stored in TGN
      instance.
      """
      return self.gpsStamp
  #End

  def createFromString(self,mString,bString,delim=" "):
      """
      This input method assumes you've opened a text file container
      and will input the data from text strings.  The assumed
      delimiter is a space but can be set in the method call.
      """
      self.delim=delim
      if self.delim == " ":
        mData=str(mString).lower().split(None)
        bData=str(bString).lower().split(None)
      else:
        mData=str(mString).lower().split(delim)
        bData=str(bString).lower().split(delim)
      mCol=mData.__len__()
      sCol=bData.__len__()
      if (mData.__len__()) != (bData.__len__()):
          raise autotrackDefineMismatchError("Array lengths %i vs %i"%(mData.__len__(),bData.__len__()))
      self.colCount=mCol
      #Break off first few elements before creating arrays
      mID=mData.pop(0)
      mDensity=mData.pop(0)
      mCount=mData.pop(0)
      mVolume=mData.pop(0)
      #
      sID=bData.pop(0)
      sDensity=bData.pop(0)
      sCount=bData.pop(0)
      sVolume=bData.pop(0)
      #
      if mID != sID:
          raise autotrackDefineIdMismatchError("Group labels do not match!")
      if mDensity != sDensity:
          raise autotrackDefineIdMismatchError("Group density values do not match!")
      if mCount != sCount:
          raise autotrackDefineIdMismatchError("Group count values do not match!")
      if mVolume!=sVolume:
          raise autotrackDefineIdMismatchError("Group volume measures do not match! %f \t %f"%(float(mVolume),float(sVolume)))
      self.__setID__(float(mID))
      self.__setDensity__(float(mDensity))
      self.__setMemberCount__(float(mCount))
      self.__setVolume__(float(mVolume))
      self.center=numpy.array([numpy.float64(j) for j in mData],'float64')
      self.bound=numpy.array([numpy.float64(j) for j in bData],'float64')
  #End createFromString

  def reverseCreateFromString(self):
    """
    Returns a tuple of two strings, that can be dumped to disk and read
    back into a TGN object with method createFromString()
    """
    stringA,stringB=self.exportToString().split("\n")
    return (stringA,stringB)
  #End reverseCreateFromString()

  def __setMemberCount__(self,mCount=0):
      """
      Sets the tally of members in the group.
      """
      self.memberCount=mCount
  #End

  def getMemberCount(self):
      """
      get the registered members listed for this TGN
      """
      return self.memberCount
  #End

  def __setVolume__(self,volume=0):
      """
      Sets the volume value of this grouping.
      """
      self.volume=volume
  #End

  def getVolume(self):
      """
      Gets the registered volume for a given TGN.
      """
      return self.volume
  #End

  def exportToString(self):
      """
      Create a text string that can be directly inserted into the
      text field of the sqlite database table TGN.
      """
      delimiter=self.delim
      outputString=""
      outputString=outputString+"%s%s"%(delimiter,self.getID())
      outputString=outputString+"%s%s"%(delimiter,self.getDensity())
      if self.colCount >= 17:
          outputString=outputString+"%s%s"%(delimiter,self.getMemberCount())
          outputString=outputString+"%s%s"%(delimiter,self.getVolume())
      for elem in self.center:
          outputString=outputString+"%s%s"%(delimiter,elem)
      outputString=outputString+"\n"
      outputString=outputString+"%s%s"%(delimiter,self.getID())
      outputString=outputString+"%s%s"%(delimiter,self.getDensity())
      if self.colCount >= 17:
          outputString=outputString+"%s%s"%(delimiter,self.getMemberCount())
          outputString=outputString+"%s%s"%(delimiter,self.getVolume())
      for elem in self.bound:
          outputString=outputString+"%s%s"%(delimiter,elem)
      return outputString
  #end exportToString()

  def getDensity(self):
      """
      Returns the value of TGN density set.
      """
      return self.density
  #End

  def __setDensity__(self,inputArg):
      """
      Sets the input density value to the TGN instance.
      """
      if type(float(0)) != type(inputArg):
          inputArg=float(inputArg)
      self.density=inputArg
  #End

  def isNULL(self):
      """
      Check the defined properties, returns true if they are all
      zeroed out which is NULL according the the matlab generator.
      """
      isNull=bool(False)
      cV=self.getCenterVector()
      bV=self.getBoundVector()
      if numpy.any(cV==0) and numpy.any(bV==0):
          isNull=bool(True)
      return isNull
  #End isNull

  def getBoundVector(self):
      """
      Gets the variance components of this TGN instance.
      """
      return self.bound
  #End getBoundVector

  def getCenterVector(self):
      """
      Gets the center of the neighborhood
      """
      return self.center
  #End getCenterVector

  def getCentroidErrorViaIndex(self,index=0):
      """
      Given an index select the that index from the Center vector
      and the bound Vector. This will return a tuple (center,bound)
      """
      centerVec=self.getCenterVector()
      boundVec=self.getBoundVector()
      vectorLength=centerVec.__len__()
      if index >= vectorLength:
          raise autotrackError("Invalid index requested exceeds\
    elements available. Elements: %i Index Sought: %i"%(vectorLength,index))
      center=centerVec[index]
      bound=boundVec[index]
      return (center,bound)
  #End getCentroidErrorViaIndex

  def isSame(self,TGN):
      """
      Checks to see if self instance is IDENTICAL to 
      TGN instance given as argument!
      """
      samePoint=bool(False)
      samePoint=self.checkOverlap(TGN,0)
      if samePoint and self.idNum==TGN.idNum and self.density==TGN.density:
          return bool(True)
      return bool(False)

  def checkOverlap(self,TGN,boundSize=0.5,mutual=bool(True)):
      """
      Check to see if SELF neighborhood overlaps with other input
      TGN class.  We define them as over lapping if they are
      withing boundSize stddevs of the center of SELF compared
      to the center of argument TGN. The last optional argument
      dictates that given the bound vector of both TGNs they should be
      able to 'walk' into each other not just A into B but that B must
      also be able to step into the neighborhood of B.
      """
      stepVectorSelf=boundSize*numpy.array(self.getBoundVector()).__abs__()
      stepVectorTGN=boundSize*numpy.array(TGN.getBoundVector()).__abs__()
      diffVector=numpy.array(
          self.getCenterVector()
          -TGN.getCenterVector()).__abs__()
      if mutual:
        if numpy.array(diffVector<=stepVectorSelf).all()\
                         and\
           numpy.array(diffVector<=stepVectorTGN).all():
          return bool(True)
        else:
          return bool(False)
      else:
        if numpy.array(diffVector<=stepVectorSelf).all():
            return bool(True)
        else:
            return bool(False)
      #End checkOverlap

  def getSeparation(self,TGN):
      """
      Gets the resultant seperation vectors and normalizes this 
      value by the boundVector, then use this to compute a 
      normalized magnitude of the vector.
      """ 
      diffVector=numpy.array(
          self.getCenterVector()
          -TGN.getCenterVector()).__abs__()
      bv=self.getBoundVector()
      sepVector=diffVector.__div__(bv)
      mag=numpy.sqrt(numpy.inner(sepVector,sepVector))
      return mag
  #End getSeperation

  def spatialAxisMinMax(self,TGN):
    """
    Return the index returned by getCenterVector for self the spatial
    axis with the smallest difference between self.getCenterVector and
    TGN.getCenterVector 
    returns a tuple of two indices {minAxisIndex,maxAxisIndex}
    """
    diffVector=numpy.array(
      self.getCenterVector()
      -TGN.getCenterVector()).__abs__()
    smallIndex=diffVector.argmin()
    largeIndex=diffVector.argmax()
    return (smallIndex,largeIndex)
  #End spatialAxisMinMax()

  def nearestTGNid(self,TGNList,boundSize=0.5):
      """
      wrapper
      """
      myN=self.nearestTGN(TGNList,boundSize)
      myID=myN.getID()
      return myID

  def nearestTGN(self,TGNList,boundSize=0.5):
      """
      Takes a list of TGNs and compares it to self
      to determine which is the closest one, this ideally
      is the same group and we can associate the group IDs. If method
      does not find TGN it returns a NONE value.
      """
      if type(list())!=type(TGNList):
          raise autotrackError("Type of input to method nearestTGNid is wrong!")
      distanceKey=list()
      for index in range(0,TGNList.__len__()):
          #Ignore entries in the list that all NULL
          if not TGNList[index].isNULL():
              if self.checkOverlap(TGNList[index],boundSize):
                  dist=self.getSeparation(TGNList[index])
                  idVal=TGNList[index].getID()
                  distanceKey.append([dist,idVal])
      distanceKey.sort()
      try:
          findID=int(distanceKey[0][1])
      except IndexError:
          findID=-1
      if findID > -1:
          for nhd in TGNList:
              if nhd.getID() == findID:
                  foundTGN=nhd
      else:
          foundTGN=None
      return foundTGN
    #End nearestTGN
#End class TGN
        
class tgnThread:
    """
    This class is a simple list object with manipulation methods which
    is called a thread that shows how the TNGs are related and it also
    tracks knots (intersections/overlaps) of two different threads.
    Thereby relating them. Three lists are tracked:
    0) The thread serial number for quick ID
    1) A thread name if given(default is thread serial number)
    2) The list of TGN what are related and sorted by time
    3) The list of knots (threads with overlaping TGNs)
    """
    def __init__(self):
      """
      HI
      """
    #End __init__()
      

class autotrackSQL:
    """
    This class provides sqlite table creation query deletion etc,
    related funtions for working the the autotrack databases.  This
    """
    def __init__(self,dbName="autotrack_default.sqlite"):
        """
        Initializes the variables associated with the autotrackSQL
        database manipulations. Setup is {TableName,Bool,{Definition}}
        """
        self.dbFile=""
        self.dbSocket=None
        #Add method to define table etc
        self.tables=dict()
        self.table['active_tgn']={'group_serial_number':'blob',
                                  'group_birthdate'    :'text',
                                  'group_gps'          :'integer',
                                  'discoverer'         :'text',
                                  'group_IFO'          :'text',
                                  'group_label'        :'text',
                                  'group_density'      :'text',
                                  'group_member_count' :'text',
                                  'group_last_seen'    :'text',
                                  'statistics'         :'blob'}

        self.table['scientist_entry']={'group_serial_number'  :'blob',
                                       'entry_date'           :'text',
                                       'scientist'            :'blob',
                                       'channels_of_interest' :'blob',
                                       'URLs_of_interest'     :'blob',
                                       'description_text'     :'text',
                                       'solution_text'        :'text',
                                       'misc_information'     :'text',
                                       'usefulness_score'     :'text',
                                       'extra_field1'         :'blob',
                                       'extra_field2'         :'blob',
                                       'extra_field3'         :'blob'}

        self.table['plot_locations']={'group_serial_number' :'blob',
                                       'line_plot'          :'blob',
                                       'aux_plots'          :'blob'}

        self.table['inactive_tgn']={'group_serial_number' :'blob',
                                      'group_birthdate'     :'text',
                                      'group_gps'          :'integer',
                                      'discoverer'          :'text',
                                      'group_label'         :'text',
                                      'group_density'       :'real',
                                      'group_member_count'  :'integer',
                                      'group_last_seen'     :'text',
                                      'statistics'          :'blob'}

        self.table['active_threads']={'thread_serial_number' :'blob',
                                      'thread_name'          :'text',
                                      'group_list'           :'blob',
                                      'knot_list'            :'blob'}

        self.table['inactive_threads']={'thread_serial_number' :'blob',
                                        'thread_name'          :'text',
                                        'group_list'           :'blob',
                                        'knot_list'            :'blob'}

        #Look for the table if it does not exist create it
        #otherwise load it into this object!
        self.__selectDB__(dbName)
        #End __init__()

    def __selectDB__(self,dbName='autotrack_default.sqlite'):
        """
        Selects the specified if it exists is reads the contents
        or if not is issues a warning and tells you the db 
        needs to be created.
        """
        #Code up how to check for a db and either create the tables
        #or read them
        #End self.__selectDB__()
        self.dbSocket=sqlite.connect(dbName)
        self.dbFile=dbName
        #Try seeing if there are tables that we want in there already
        try:
            self.dbSocket.execute("select * from %s"%(self.defineTables[0][0]))
        except:
            sys.stdout.write("It appears that we need to create the tables.\n")
            self.createTables()
        #End __selectDB__()

    def DEPRICATEDgetTableIndex(self,name=""):
        """
        Given a string searches the tuple self.defineTables to
        determine the index of that table to various function calls.
        """
        answer=int(-1)
        for index in range(self.defineTables.__len__()):
            if self.defineTables[index][0].lower() == name.lower():
                answer=index
        if not(type(int()) == type(answer)):
            raise autotrackError("getTableIndex type invalid!\n %s"(answer))
        return answer
    #End getTableIndex

    def __DEPRICATEDcreateSingleTable__(self,name):
        """
        This method will generate a table from the list of tables
        and definitions specified by self.autotrackTableDef or
        self.auxTableDef
        """
        tableIndex=self.DEPRICATEDgetTableIndex(name)
        thisTable=self.defineTables[tableIndex]
        tableName=thisTable[0]
        required=thisTable[1]
        colDefs=thisTable[2]
        rowString=""
        for col in colDefs:
            rowString=rowString+"%s "%(col)
        commandString="create table %s (%s)"%(tableName,rowString)
        try:
            self.dbSocket.execute(commandString)
        except:
            sys.stderr.write("Error trying to create table.\n")
            self.dbSocket.commit()
            return
        self.dbSocket.commit()
        #End __createSingleTable__()

    def __createSingleTable__(self,name=None,tabledef=None):
        """
        This method will generate a table from the list of tables
        and definitions specified by self.autotrackTableDef or
        self.auxTableDef
        """
        if (name == None or tabledef == None):
            return
        tableName=name
        colDefs=tabledef
        rowString=""
        for col in colDefs:
            rowString=rowString+"%s "%(col)
        commandString="create table %s (%s)"%(tableName,rowString)
        try:
            self.dbSocket.execute(commandString)
        except:
            sys.stderr.write("Error trying to create table.\n")
            self.dbSocket.commit()
            return
        self.dbSocket.commit()
        #End __createSingleTable__()
        
    def getSocket(self):
        """
        Returns an object which is a raw handle to the dqlite DB.
        This is equivalent to being returned X from
        X=sqlite.connect('dbfilename')
        """
        return self.dbSocket
        #End getSocket

    def createTables(self,overwrite=bool(False)):
        """
        This function call will create the specified sqlite db unless
        there exists one already.  In that case we throw an error
        unless we explicitly want to overwrite that table.
        """
        for index in self.tables.keys():
            sys.stdout.write("Creating table %s\n"%(index))
            self.__createSingleTable__(index,self.tables[index])
        self.dbSocket.commit()
        #End createTables()

    def getThreads(self,table='tgn',time=None):
        """
        This method will fetch the groups defined at time NOW in table
        TGN(default).  It returns these groups as a list of TGN objects.
        """
        #Search the group threads and grab a list of all threads.  
        #Take head of each thread and try for matching it against newly
        #seen TNGs to increase the thread length.  
        #
        print "HI"
        #End getThreads()

    def getThreadHeads(self,table='tgn',time=None):
        """
        This returns only the head of the thread to keep memory
        requirements to a minimum.  The thread object will only
        be loaded if a manipulation to the thread or thread report
        method is called.
        """
        print "Getting all TGN heads from threads."
        #End getThreadHeads()

    def getThreadByName(self,table='tgn',name=None):
        """
        This function will return a thread object (list) for a 
        TGN thread that exists in the table specified.  By default
        we search the active threads in table TGN.
        """
        print "Getting single thread by name"
        #End 

    def getThreadHeadByName(self,table='tgn',name=None):
        """
        Returns a single TGN object the thread head for quoted
        thread.
        """
        print "Getting TGN head from thread"
        #End getThreadHeadByName

    #def loadThread(self,serialNumber=None)
    #def syncThread(self)
    #def deactivateThread(self)
    #def activateThread(self)
    #def deleteThread(self)

    
