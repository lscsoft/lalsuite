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
__date__ = '$Date$'
__version__ = ''


import os
from numpy import ndarray,dtype,NaN,Inf,zeros,mean,sin,\
     pi,all,any,array,arange,inner,cross,isnan,arcsin,\
     mod
from pylal.metaarray import TimeSeries,TimeSeriesList
import copy
import sys
from StringIO import StringIO
from commands import getstatusoutput
from scipy import interp
from followup_utils import connectToFrameData
from pylal import metaarray
import warnings

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
This file will hold all the methods required to run the  c^2t(Cristina
Cesar Triangulation) DetChar tool.
"""


#
# Iterate in matrix computing row dot vectors
#
def CSTdot(A,B):
    """
    Input two matrices Nx3 in size.
    Calculates the inner product row by row.
    """
    if( A.shape != B.shape) and A.ndim > 2:
      os.abort()
    if A.ndim == 2:
      return array([inner(A[i],B[i]) for i in arange(A.shape[0])])      
    else:
      return array(inner(A,B))
#
# Iterate in matrix computing row cross vectors
#
def CSTcross(A,B):
    """
    Input two matrices Nx3 in size.
    Calculates the cross product row by row.
    """
    if( A.shape != B.shape) and A.ndim > 2:
      os.abort()
    if A.ndim == 2:
      return array([cross(A[i],B[i]) for i in arange(A.shape[0])])
    else:
      return array(cross(A,B))

#
# Intersect two planes
#
def intersectPlanes(normA,normB,coordA,coordB,wiggleA=0.0):
    """
    Expects a 3xM matrices of 3D vectors for the plane norms
    also in addition to that we take vectors from the "origin"
    (0,0,0) (single 3x1 matrix) to be the coordinate input.  The
    output is a tuple of 3xM arrays parameterized by their index,
    intersection vector and point on that line. Output is
    (linePoint,LineVector) ((3xM),(3xM)). The input wiggleA is the
    wiggle angle in radian to consider two vectors parrallel.
    """
    #If inputs a 'lone' vectors make them a 3x1 matrix!
    #Normalize the input vectors, but NOT coords!
    for j in arange(len(normA)):
        normA[j]=normA[j]/CSTdot(normA[j],normA[j])
        normB[j]=normB[j]/CSTdot(normB[j],normB[j])
    #Cycle through all the rows of input matrix
    lineVec=zeros(normA.shape)
    linePoi=zeros(normA.shape)
    for i in arange(len(normA)):
        #Determine if planes parallel?
        if arcsin(CSTdot(normA[i],normB[i])) >= (1-wiggleA):
            #Set values to NaN for given 'i'
            lineVec[i]=NaN
            linePoi[i]=NaN
        else:
            #Find normal vector along line of intersection!
            tmpVec=CSTcross(normA[i],normB[i])
            lineVec[i]=tmpVec/CSTdot(tmpVec,tmpVec)
            print "Intersecting line vector(unit) %s"%(lineVec)
            #Find point along line of intersection!
            if isnan(lineVec[i]).any() or all(lineVec[i]==0.0):
              print "Vector of plane intersection line is ZERO Vector!"
              lineVec[i]=NaN
              linePoi[i]=NaN
            else:
              #Determine which zero axis vector points through!
              zeroIndex=lineVec.argmax()
              #Map out x,y,z to index of i,j,k with lowest index as 'i'
              (iX,iY,iZ)=mod(range(zeroIndex,zeroIndex+3),3).tolist()
              print "Vector elements:%s,%s,%s"%(iX,iY,iZ)
              #(x,y,z)==(l0,l1,l2)
              linePoi[i][iZ]=0
              linePoi[i][iX]=(CSTdot(normB[i],coordB)-(CSTdot(normA[i],coordA)*normB[i][iY])/normA[i][iY])/(normB[i][iX]-(normB[i][iY]*normA[i][iX])/normA[i][iY])
              linePoi[i][iY]=(CSTdot(normA[i],coordA)-normA[i][iX]*linePoi[i][iX])/normA[i][iY]
    return (lineVec,linePoi)

#
# Intersect two lines
#
def intersectLines(vecA,vecB,coordA,coordB,wiggleA=0.0):
    """
    Expects inputs of 3xM matrices representing vectors along a line
    and a vector 3xM of matching points on lines for those vectors.
    Output is in a 3xM collection of intersecting points.
    Examples:
    Input
    ([[vX,vY,vZ]],[[lX,lY,lZ]])

    Output
    [[i,j,k]]
    A set of lines is considered parallel if angle(rads) between
    then is less than wiggleA, in this case the output coordinate
    is (NaN,NaN,NaN).
    """
    #Error check input shapes etc!
    #Normalize the input vectors, but NOT coords!
    for j in arange(len(vecA)):
        vecA[j]=vecA[j]/CSTdot(vecA[j],vecA[j])
        vecB[j]=vecB[j]/CSTdot(vecB[j],vecB[j])
    #Generate a matrix of triangulated coordinates
    triPoi=zeros(vecA.shape)
    for i in arange(len(vecA)):
        #Check if line vectors are parallel?(vecA&vecB unit vectors)
        if arcsin(CSTdot(vecA[i],vecB[i]))<=wiggleA:
          #Set values to NaN for given 'i'
          triPoi[i]=NaN
        else:
          (iX,iY,iZ)=(0,1,2)
          #Find intersection of the two lines.
          #P(s)=P_0+s*\vec{P}
          gammaConst=(\
            coordA[i][iY]*vecA[i][iZ] - vecA[i][iZ]*coordB[i][iY] +\
              vecA[i][iY]*(coordB[i][iZ]-coordA[i][iZ])\
            )/\
            (vecA[i][iY]*vecB[i][iZ]+vecA[i][iZ]*vecB[i][iY])
          triPoi=coordB+(gammaConst*vecB)
    return triPoi

class cSqrTSensor(object):
  """
  This class is just method globber that works expects as input
  TimeSeries objects defined by pylal.frutils.  This class only
  works in a cartesian basis!  The distance should be measure in units
  of meters.
  """
  def __init__(self,\
               name="Unknown",\
               xAxis=None,\
               yAxis=None,\
               zAxis=None,\
               coordinate=(None,None,None),\
               swapSign=(False,False,False),\
               tType="normal"\
               ):
    """
    This object is initialized with 3 data streams (1xM) which should be the
    intesities of each spatial coordinate seen by the sensor in
    question.  In addition to this information, the cartesian
    coordinates relative to some origin should be specified.  Inside
    the class the data will be represented as collections of row
    vectors in a 3xM matrix.  The units
    are not checked between different cSqrTSensor objects.  The user
    should take care to make sure this is sensible!  
    """
    self.name=name.strip().upper()
    self.__tType__=tType.strip().lower()
    if self.__tType__ != "normal":
      os.abort()
    #
    # Set variable type as REAL8
    self.dtype='float64'
    #
    # Set angular wiggle room to define parralel
    wiggleAngle=3.0 #Degrees
    wiggleDistance=25.0 #Meters
    self.__wiggleA__=None
    self.__wiggleD__=None
    self.setParallelMismatchAngle(wiggleAngle)
    self.setHypersphereMismatch(wiggleDistance)
    #
    # Split off the metadata if input is a TimeSeries
    #
    if isinstance(xAxis,TimeSeries):
      self.metaX=xAxis.metadata
    else:
      self.metaX=None
    if isinstance(yAxis,TimeSeries):
      self.metaY=yAxis.metadata
    else:
      self.metaY=None
    if isinstance(zAxis,TimeSeries):
      self.metaZ=zAxis.metadata
    else:
      self.metaZ=None
    #
    # Verify Input data properties if TimeSeries
    #
    if all(self.metaX,self.metaY,self.metaZ):
      if not (
        xMetaDict["segments"].intersects_all(yMetaDict["segments"]) \
        and \
        yMetaDict["segments"].intersects_all(zMetaDict["segments"]) \
        and \
        zMetaDict["segments"].intersects_all(xMetaDict["segments"]) \
        ):
        raise Exception, \
              "Inputs for X,Y,Z data don't cover identical segments!"
      #
      # Check sampling rates! Add ability to resample on the fly!
      #
      if not(\
        xMetaDict["dt"] == \
        yMetaDict["dt"] == \
        zMetaDict["dt"] \
        ):
        raise Exception, \
              "Input data for sampling rates on data is not consistent!"
    else:
      sys.stdout.write("Inputs are not TimeSeries ojects can not verify data traits.\n")

    #
    # Save the coordinate information
    # (Coordinate 3x1 Array [X,Y,Z] NOT [[X,Y,Z]])
    if isinstance(coordinate,ndarray):
      self.coordinate=coordinate
    else:
      self.coordinate=array(coordinate,dtype=self.dtype)
    #
    # Save the input data
    #
    self.data=zeros((len(xAxis),3),self.dtype)
    self.data[:,0]=xAxis
    self.data[:,1]=yAxis
    self.data[:,2]=zAxis
    #
    # End Init
    #

  def setParallelMismatchAngle(self,\
                               wiggleAngle=3.0\
                               ):
    """
    Enter an angle in degrees to define a tolerance on the
    parallelness of two lines for determining potential
    intersection of lines and planes.
    """
    self.__wiggleA__=sin(wiggleAngle/(2.0*pi))#Radians
    
  def setHypersphereMismatch(self,\
                             wiggleDistance=25.0\
                             ):
    """
    Enter an angle in degrees to define a tolerance on the
    parallelness of two lines for determining potential
    intersection of lines and planes.
    """
    self.__wiggleD__=wiggleDistance

  def triangulate(self,\
                  otherSensor=[None],\
                  ):
    """
    Wrapper to pick the correct triangulation machinery.
    """
    if not isinstance(otherSensor,list):
      otherSensor=[otherSensor]
    if not all([isinstance(x,cSqrTSensor) for x in otherSensor]):
      raise Exception, "Expected list of cSqrTSensor objects, got list of other objects!"
    if self.__tType__ == "normal":
      return self.__normalTriangulate__(otherSensor)
    else:
      os.abort()

  def __normalTriangulate__(self,\
                            otherSensor=[None],\
                            ):
    """
    Assuming the input is normals to 3D planes in 3 space we will
    intercept the planes resulting in either points in 3D space or 2D
    lines whichever is the most physically interesting.
    """
    #
    # Check for name space collisions
    if self.name in [x.name for x in otherSensor]:
      raise Exception, "Namespace collision between sensors to triangulate!"
    #
    # Perform triangulation
    # (self,A,B)
    planeIntersections=dict()
    pKeyMask="%s_%s"
    if len(otherSensor) >= 1:
      #Intersect planes pairwise
      sensorList=list()
      sensorList.append(self)
      sensorList.extend(otherSensor)
      #Identify uniq pairs
      for aa,A in enumerate(sensorList):
        for bb,B in enumerate(sensorList):
          myKey=pKeyMask%(A.name,B.name)
          symKey=pKeyMask%(B.name,A.name)
          if (aa!=bb) and \
             not (myKey in planeIntersections.keys()) and \
             not (symKey in planeIntersections.keys()):
            planeIntersections[pKeyMask%(A.name,B.name)]=intersectPlanes(A.data,\
                                                                         B.data,\
                                                                         A.coordinate,\
                                                                         B.coordinate,
                                                                         wiggleA=self.__wiggleA__)
    else:
      raise Exception, "Need at least two normals to intersect planes."
    #In case of only on pair of plane given as inputs
    if len(planeIntersections.keys()) == 1:
      myKey=planeIntersections.keys()
      return planeIntersections[myKey[0]]
    #If we have a pairs of lines to intersect
    lineIntersections=dict()
    for (nameA,lineA) in planeIntersections.iteritems():
      for (nameB,lineB) in planeIntersections.iteritems():
        lineIntersections[pKeyMask%(nameA,nameB)]=intersectLines(lineA[0],\
                                                                 lineB[0],\
                                                                 lineA[1],\
                                                                 lineB[1],\
                                                                 wiggleA=self.__wiggleA__)
    #
    # If points at time index all inside this hypersphere return mean point!
    #
    resultPoint=zeros(lineIntersections[0],dtype=self.dtype)
    resultVectors=zeros(lineIntersections[0],dtype=self.dtype)
    #If points close enough return centroid else None
    for i in arange(len(lineIntersections[0])):
      tmpCentroid=zeros(lineIntersections[0],dtype=self.dtype)
      for j in arange(len(lineIntersections)):
        tmpCentroid+=lineIntersections[j]
      tmpCentroid=tmpCentroid/float(len(lineIntersections))
      resultPoint[i]=tmpCentroid
      if max([sqrt(CSTdot(x-tmpCentroid,x-tmpCentroid)) \
              for x in lineIntersections.iteritems]) > \
              self.__wiggleD__:
        resultPoint[i]=NaN
      else:
        resultPoint[i]=tmpCentroid
    #
    # Return triangulation coordinates
    #
    return (resultVectors,resultPoint)
