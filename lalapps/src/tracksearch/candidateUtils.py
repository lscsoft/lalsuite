#!/usr/bin/env python2.3

__author__ = 'Charlie Torres <charlie@phys.utb.edu>'
__date__ = '$Date$'
__version__ = ''

import getopt
import math
import os
import string
import sys
import time
import copy


"""
This file provides the class and methods needed to process candidate
files. These claases will provide three features.
Globbing -> Condense list of candidates down into a single larger file.
Clobber -> take two files with dissimiliar or similiar size maps
           processed removing any curves from A that live in list B.
Snapshot -> create a PGM image to use a mental snapshot of where the
            curves found reside in the TF space
"""
class kurve:
    """
    This is a class that will hold the data in primitives of python and
    provide higher level methods to questioning the curve object.
    """
    def __init__(self,curveId,length,power):
        self.curveId=int(curveId)
        self.length=int(length)
        self.power=float(power)
        self.element=list()
    #End __init__ method

    def __len__(self):
        """
        Determine pixels in curve as an integer and returns that number.
        """
        return self.element.__len__()
    #End __len__ method
    
    def appendPixel(self,Row,Col,gpsStamp,Freq,Power):
        """
        Add a new pixel to our curve.  This should not be used manually in
        most cases.  It is invoked via the loadfile method in the candidateList
        class
        """
        self.element.append([Row,Col,gpsStamp,Freq,Power])
        self.length=self.__len__()
    #End method appendPixel

    def printCurveID(self):
        """
        Method to return the curveID of the structure IDs are not required
        to be unique to the list.
        """
        return self.curveId
    #End printCurveID method
    
    def printStartGPS(self):
        """
        This method will return the start GPS time as a text string
        """
        if self.element.__len__() > 0:
            result=self.element[0][2].display()
        else:
            print "Object appears to be empty!"
            result=gpsInt(0,0)
        return result
    #End method print StartGPS

    def startGPS(self):
        """
        This method will return the start GPS time as a text string
        """
        if self.element.__len__() > 0:
            result=self.element[0][2]
        else:
            print "Object appears to be empty!"
            result=gpsInt(0,0)
        return result
    #End method startGPS

    def printStopGPS(self):
        """
        This method will return the start GPS time as a text string
        """
        curveSize=self.element.__len__()
        if curveSize > 0:
            result=self.element[curveSize-1][2].display()
        else:
            print "Object appears to be empty!"
            result=gpsInt(0,0)
        return result
    #End method print StopGPS

    def stopGPS(self):
        """
        This method will return the start GPS time as a text string
        """
        curveSize=self.element.__len__()
        if curveSize > 0:
            result=self.element[curveSize-1][2]
        else:
            print "Object appears to be empty!"
            result=gpsInt(0,0)
        return result
    #End method stopGPS

    def printStartFreq(self):
        """
        This method will print the start Frequency of the curve
        """
        curveSize=self.element.__len__()
        if curveSize > 0:
            result=self.element[0][3]
        else:
            print "Object appears to be empty!"
            result=float(-1)
        return result
    #End method printStartFreq

    def printStopFreq(self):
        """
        This method will print the stop Frequency of the curve
        """
        curveSize=self.element.__len__()
        if curveSize > 0:
            result=self.element[curveSize-1][3]
        else:
            print "Object appears to be empty!"
            result=float(-1)
        return result
    #End method printStopFreq
        
    def getCoordTriple(self,n):
        """
        This method will get the GPStime,Freq and Power of a pixel with
        index value n and return it as a three element tuple
        """
        curveSize=self.element.__len__()
        if ((n >= curveSize) and (n < 0)):
            print "Requested Index exceeds available data! :",n,curveSize
            return gpsInt(-1,0),float(-1),float(-1)
        else:
            return self.element[n][2],self.element[n][3],self.element[n][4]
    #End method getCoordTriple

    def getCoordBins(self,n):
        """
        This method will get the interger bin numbers associated with the
        curve.  These values are rather useless once your globbed together
        many coordinateList structures.
        """
        curveSize=element.__len__()
        if ((n > curveSize) and (n >= 0)):
            print "Requested Index exceeds available data! :",n
            return int(-1),int(-1)
        else:
            return element[n][0],element[n][1]
    #End method getCoordBins

    def getKurveHeader(self):
        """
        This method will return the values as n primatives to the calling
        routine: curveId,length,power
        """
        return self.curveId,self.length,self.power
        #End method getKurveHeader

    def getKurveDataBlock(self):
        """
        This method will return an copied iterable N*5 2 dimensional list.
        This method is normally invoked by writefile method in
        candidateList class
        """
        return copy.deepcopy(self.element)
        #End method getKurveDataBlock
        
class gpsInt:
    """
    Provides method for converting our string format to GPS of X,Y to
    a long int with implicit decimal for basic numerical manipulations
    """

    def __init__(self,gpsSeconds,gpsNanoSeconds):
        """
        Setup an instance of gpsfloat for use in the code
        """
        try:
            self.gpsSeconds=int(gpsSeconds)
            self.gpsNanoSeconds=int(gpsNanoSeconds)
        except ValueError:
            foundType1=type(gpsSeconds)
            foundType2=type(gpsNanoSeconds)
            print "Error can not cast input arguments to integer types."
            print "Found types :",foundType1,foundType2
            print "Values      :",gpsSeconds,gpsNanoSeconds
            os.abort()
    #End init method

    def __makeInt__(self):
        secPartIn=str(self.gpsSeconds)
        nanoPartIn=str(str(self.gpsNanoSeconds).rjust(9)).replace(' ','0')
        intIn=int(secPartIn+nanoPartIn)
  #        print "As Ints: ",secPartIn,nanoPartIn,intIn
        return intIn
    #End __makeInt__ method

    def __splitInt__(self,inputInt):
        negFlag=False
        if inputInt < 0:
            inputInt = abs(inputInt)
            negFlag=True
        strInput=str(inputInt).zfill(18)
        gpsSecondInput=int(strInput.__getslice__(0,9))
        gpsNanoInput=int(strInput.__getslice__(9,18))
        if negFlag == True:
            gpsSecondInput=-1*gpsSecondInput
        return gpsSecondInput,gpsNanoInput
    #End __splitInt__ method
    
    def __add__(self,inTime):
        """
        Add the argument gpsInt class to the calling gpsInt object ie
        x.__add__(y) like the operation x+y
        """
        inTimeInt=inTime.__makeInt__()
        selfTimeInt=self.__makeInt__()
        intAnswer=inTimeInt+selfTimeInt
        gpsSecondAnswer,gpsNanoAnswer=self.__splitInt__(intAnswer)
        return gpsInt(gpsSecondAnswer,gpsNanoAnswer)
    #End __add__ method

    def __sub__(self,inTime):
        inTimeInt=inTime.__makeInt__()
        selfTimeInt=self.__makeInt__()
        intAnswer=selfTimeInt.__sub__(inTimeInt)
        gpsSecondAnswer,gpsNanoAnswer=self.__splitInt__(intAnswer)
        return gpsInt(gpsSecondAnswer,gpsNanoAnswer)
    #End __sub__ method

    def __div__(self,divisor):
        secPart=self.gpsSeconds.__floordiv__(divisor)
        secRemain=self.gpsSeconds.__mod__(divisor)
        nanoSecPart=int(self.gpsNanoSeconds.__floordiv__(divisor)).__add__(secRemain)
        return gpsInt(secPart,nanoSecPart)
    #End __div__ method
    
    def display(self):
        secPartIn=str(self.gpsSeconds)
        nanoPartIn=str(str(self.gpsNanoSeconds).rjust(9)).replace(' ','0')
        result=secPartIn+'.'+nanoPartIn
        return result
    #End display method
    
    def __diskPrint__(self):
        secPartIn=str(self.gpsSeconds)
        nanoPartIn=str(str(self.gpsNanoSeconds).rjust(9)).replace(' ','0')
        result=secPartIn+','+nanoPartIn
        return result
<<<<<<< candidateUtils.py
    #End __diskPring__ method

    def __abs__(self):
        if (self.gpsSeconds < 0) and (self.gpsNanoSeconds < 0):
            print "Major weirdness with data! See ",self.gpsSeconds,self,gpsNanoSeconds
            print "If you see this something went really wrong!"
            os.abort()
        return gpsInt(abs(self.gpsSeconds),abs(self.gpsNanoSeconds))
=======
    #End __diskPring__ method
>>>>>>> 1.5
#End gpsInt class

class candidateList:
    """
    Provides basic IO for manipulation of candidate lists
    """

    def __init__(self):
        self.totalCount=int(0)
        self.filename=list()
        self.numFbins=int(-1)
        self.numTbins=int(-1)
        self.gpsSpan=gpsInt(-1,0)
        self.freqSpan=float(-1)
        self.freqWidth=float(0)
        self.gpsWidth=gpsInt(0,0)
        self.sorted=False
        #Should be list of objects of class kurve
        self.curves=[]
    #End init method

    def __sec2pixel__(self,seconds):
        """
        This method returns the number of rounded pixels which represents
        the requested number of seconds given via:
        numOfPixel=ceil(seconds/pixelDuration)
        """
        return int(math.ceil(seconds/self.gpsWidth.__asFloat__()))

    def __filemaskGlob__(self):
        """
        Method to use the filename is self.filename list to create a new glob
        file name to write to the disk with.
        """
        if self.filename.__len__() == 0:
            print 'Error with header information in instance'
            return 'DefaultName.file'
        else:
            filename='GLOB:'+os.path.basename(self.filename[0])+\
            'Entries:'+str(self.filename.__len__())
        return filename

    def __filemaskClob__(self,clobWith):
        """
        Method to us if we want a automagic filename to write the clobber
        result file to.
        """
        clobberVictim=str(str(os.path.basename(self.filename[0])).split('.')[0])
        clobberPerp=str(os.path.basename(clobWith.filename[0])).split('.')
        filename='CLOBBER:'+clobberVictim+clobberPerp

    def loadfile(self,inputFilename):
        """
        Reads in a candidate list from the disk with a given filename
        """
        input_fp=open(inputFilename,'r')
        content=input_fp.readlines()
        input_fp.close()
        if content.__len__() < 1:
            print "Error no lines in file?"
            print "Check :",inputFilename
            os.abort()
        self.totalCount=int(list(content[0].split(':'))[1])
        self.filename=[inputFilename]
        #Check for comment lines -> #Comment
        index2pop=[]
        for index in range(0,int(content.__len__())):
            if content[index].startswith('#'):
                index2pop.append(index)
        index2pop.reverse()
        for index in index2pop:
            content.pop(index)
        #Check for empty lines
        index2pop=[]
        for index in range(0,int(content.__len__())):
            if content[index].startswith('\n'):
                index2pop.append(index)
        index2pop.reverse()
        for index in index2pop:
            content.pop(index)
        #Simple file parity check
        if ((int(content.__len__()).__mod__(2) != 0)):
            print 'Error parsing :',self.filename," Number of data lines :",content.__len__()
            os.abort()
<<<<<<< candidateUtils.py
        if (content.__len__() > 1) and (self.totalCount == 0):
            print "Hum?  The file appears inconsistent."
            print "File has ",content.__len__(),"lines with header listing ",self.totalCount," entries."
            print inputFilename
=======
        if (content.__len__ > 0) and (self.totalCount == 0):
            print "Hum?  The file appears inconsistent."
            print inputFilename
>>>>>>> 1.5
        if (content.__len__ > 0):
             walkIndex=0
             while walkIndex < content.__len__():
                 entry=[]
                 entryHeader=[]
                 entryElement=[]
                 entryElementTXT=[]
                 entryHeader=str(list(content[walkIndex].split(':'))[1]).split(',')
                 self.curves.append(kurve(entryHeader[0],entryHeader[1],\
                                          entryHeader[2]))
                 walkIndex=walkIndex.__add__(1)
                 entryElementTXT=str(content[walkIndex].replace(';',',')).split(':')
                 walkIndex=walkIndex.__add__(1)
                 entryList=[]
                 #Parse the text line with our delimiters which designate
                 #pixels in the curve appending them to the curve object
                 for coord in entryElementTXT:
                     tmpElement=str(coord).split(',')
                     self.curves[self.curves.__len__()-1].appendPixel(\
                         int(tmpElement[0]),int(tmpElement[1]),\
                         gpsInt(tmpElement[2],tmpElement[3]),\
                         float(tmpElement[4]),float(tmpElement[5]))
<<<<<<< candidateUtils.py
             #Determine the bin widths in this structure
             if self.totalCount > 0:
                 self.findBinWidths()
             if self.totalCount != self.curves.__len__():
                 print "Possible problem, Inconsistent file :",inputFilename
=======
                 #Determine the bin widths in this structure
                 if self.totalCount > 0:
                     self.findBinWidths()
>>>>>>> 1.5
        else:
            print "No candidate entries found in:",inputFilename
    #End loadfile method

    def writefile(self,outputFilename):
        """
        Write the current candidate list structure to disk using given
        filename
        """
        output_fp=open(outputFilename,'w')
        textLine="#Total Curves:"+str(self.totalCount)+"\n"
        output_fp.write(textLine)
        output_fp.write('# Legend: Col,Row;gpsSec,gpsNanoSec,Freq,depth\n')
        output_fp.write('#\n')
        for entry in self.curves:
            CurveId,Length,Power=entry.getKurveHeader()
            text="Curve number,length,power:"+str(CurveId)+','+str(Length)+','+str(Power)+'\n'
            output_fp.write(text)
            text=""
            for elements in entry.getKurveDataBlock():
                text=text+str(elements[0])+','+str(elements[1])+\
                      ';'+str(elements[2].__diskPrint__())+','\
                      +str(elements[3])+','+str(elements[4])+':'
            text=text.rstrip(':')+'\n'
            output_fp.write(text)
    #End writefile method

    def sortList(self):
        """
        This sorts the input curves into an order with lowest gpsStart time
        first.  The sorting keys only on the curves start time.
        """
        #Created keye list gpsInt.display(),OriginalIndex
        currentIndex=0
        keyedList=[]
        for entry in self.curves:
            tmpPair=[]
            tmpPair.append((entry.startGPS()).__makeInt__())
            tmpPair.append(currentIndex)
            currentIndex=currentIndex.__add__(1)
            keyedList.append(tmpPair)
        #Sort this keyed list lowest->highest
        keyedList.sort()
        newList=[]
        for entry in keyedList:
            newList.append(self.curves[entry[1]])
        if newList.__len__() != self.curves.__len__():
            print "Error with sorting routine!"
            os.abort()
        #Reassign the list to newly sorted version
        self.sorted=True
        self.curves=copy.deepcopy(newList)
    #End sortList method

    def findBinWidths(self):
        """
        Use the information loaded to guess the bin size T and F
        we will not know the original map bandwidth but we
        won't really need it Col->F Row->T
        """
        TDArrayFreq=[]
        TDArrayGPS=[]
        for curveElement in self.curves:
            for entry in curveElement.getKurveDataBlock():
                TDArrayFreq.append([entry[0],entry[3]])
                TDArrayGPS.append([entry[1],entry[2]])
#        print "We will use ",TDArrayFreq.__len__()," total pixels to determine time bin width and freq bin width"
        TDArrayFreq.sort()
        TDArrayGPS.sort()
        freqSum=float(0)
        freqSumCount=0
        for i in range(1,TDArrayFreq.__len__()):
            if int(TDArrayFreq[i][0]) - int(TDArrayFreq[i-1][0]) == 1:
                AddMe=float(TDArrayFreq[i][1])-float(TDArrayFreq[i-1][1])
                AddMe=abs(AddMe)
                freqSum=freqSum.__add__(AddMe)
                freqSumCount=freqSumCount+1
        if freqSumCount != 0:
            avgFreqWidth=freqSum.__div__(freqSumCount)
        else:
            avgFreqWidth=0
        gpsSum=gpsInt(0,0)
        gpsSumCount=0
        for i in range(1,TDArrayGPS.__len__()):
            if int(TDArrayGPS[i][0]) - int(TDArrayGPS[i-1][0]) == 1:
                AddMe=TDArrayGPS[i][1].__sub__(TDArrayGPS[i-1][1])
                AddMe=AddMe.__abs__()
#                print TDArrayGPS[i][1].display(),TDArrayGPS[i-1][1].display(),AddMe.display()," As INTS ",TDArrayGPS[i][1].__makeInt__(),TDArrayGPS[i-1][1].__makeInt__(),AddMe.__makeInt__()
#                print "Showing  i,AddMe ",i,AddMe.gpsSeconds,AddMe.gpsNanoSeconds,AddMe.display()
#                print "Index ",i," of ",TDArrayGPS.__len__()," Adding value : ",AddMe.display()," To ",gpsSum.display()
                gpsSum=gpsSum.__add__(AddMe)
                gpsSumCount=gpsSumCount+1
<<<<<<< candidateUtils.py
        if gpsSumCount > 0:
=======
        print "Sum Count: ",gpsSumCount,"Total GPS Time: ",gpsSum.display()
        if gpsSumCount > 0:
>>>>>>> 1.5
            avgGpsWidth=gpsSum.__div__(gpsSumCount)
        else:
            avgGpsWidth=gpsInt(0,0)

        self.freqWidth=float(avgFreqWidth)
        self.gpsWidth=avgGpsWidth
    #End findBinWidths method

    def globList(self,inputCandidateList):
        """
        Take list arguement and concatinate it with the self candidate list
        only if the frequency and time bin widths match
        """
        iCL=inputCandidateList
        globList=copy.deepcopy(self)
        globTolerance=1e-5
        gpsDiff=self.gpsWidth.__sub__(iCL.gpsWidth)
        if ((float((self.freqWidth-iCL.freqWidth)).__abs__() < globTolerance)\
           and \
           (float(gpsDiff.display()).__abs__() < globTolerance)):
            globList.totalCount=self.totalCount.__add__(iCL.totalCount)
            globList.filename.extend(iCL.filename)
            globList.curves.extend(iCL.curves)
            return globList
        elif ((self.freqWidth==0) or (int(self.gpsWidth.__makeInt__())==0)):
            globList.freqWidth=self.freqWidth
            globList.gpsWidth=self.gpsWidth
            globList.totalCount=self.totalCount.__add__(iCL.totalCount)
            globList.filename.extend(iCL.filename)
            globList.curves.extend(iCL.curves)
            return globList
        elif ((iCL.freqWidth==0) or (int(iCL.gpsWidth.__makeInt__())==0)):
            globList.freqWidth=iCL.freqWidth
            globList.gpsWidth=iCL.gpsWidth
            globList.totalCount=self.totalCount.__add__(iCL.totalCount)
            globList.filename.extend(iCL.filename)
            globList.curves.extend(iCL.curves)
            return globList
        else :
            print "Can not glob lists, due to inconsistent values!"
            print self.freqWidth,'VS',iCL.freqWidth
            print self.gpsWidth.display(),'VS',iCL.gpsWidth.display()
            print "Returning empy list"
            return candidateList()
    #End globList method

    def clusterClobberWith(self,inputReferenceList):
        """
        Take instance and compare candidate entries against input reference.
        Delete from this instance any entries which fit into candidate entries
        located in the inputReferenceList given to this method.  Return the
        a reduced candidate reference instance.
        """
        iRL=copy.deepcopy(inputReferenceList)
        #Check to see if list has been sorted if not do it.
        if (self.sorted != True):
            self.sortList()
        if (iRL.sorted != True):
            iRL.sortList()
        selfStopIndex=self.totalCount
        iRLStopIndex=iRL.totalCount
        #Do the two candidate lists overlap in time?
        selfStart=(self.curves[0].startGPS()).__makeInt__()
        selfStop=(self.curves[self.curves.__len__()-1].stopGPS()).__makeInt__()
        iRLStart=(self.curves[0].startGPS()).__makeInt__()
        iRLStop=(self.curves[self.curves.__len__()-1].stopGPS()).__makeInt__()
        if not((selfStop >= iRLStart) or (iRLStop >= selfStart)):
            print "It appears as if these two candidateList instances do not have overlapping data!"
            print "Self's    Boundaries: ",selfStart,selfStop
            print "Clobber's Boundaries: ",iRLStart,iRLStop
            print "Returning copy of the original instance."
            return copy.deepcopy(self)
        #Since lists overlap start processing candidates.
        #Setup bin sizes
        hgw=halfGPSwidth=iRL.gpsWidth.__div__(2)
        hfw=halfFreqWidth=iRL.freqWidth.__div__(2)
        candidateToRemove=[]
        trialCount=0
        for i in range(0,self.totalCount):
            #from kurve i find all j curves that equal or span that time interval
            #build a list of those curves indexs to analyze
            #print "Checking i index:",i
            iRLmatchIndex=[]
            for j in range(0,iRL.totalCount):
                if (iRL.curves[j].startGPS().__makeInt__() <= self.curves[i].startGPS().__makeInt__())\
                   and \
                   (iRL.curves[j].stopGPS().__makeInt__() >= self.curves[i].stopGPS().__makeInt__()):
                    iRLmatchIndex.append(j)
            for j in iRLmatchIndex:
                #print "Curve Indexes self,iRL",self.curves[i].printCurveID(),iRL.curves[j].printCurveID(),"Index :",i,j
                #Build the self TF list
                selfCurve=[]
                for index in range(0,self.curves[i].__len__()):
                    tmpT,tmpF,JNK=self.curves[i].getCoordTriple(index)
                    selfCurve.append([tmpT,tmpF])
                #Build the iRL TF list
                irlCurve=[]
                for index in range(0,iRL.curves[j].__len__()):
                    tmpT,tmpF,JNK=iRL.curves[j].getCoordTriple(index)
                    irlCurve.append([tmpT,tmpF])
                boundPixels=0
                kStop=selfCurve.__len__()
                lStop=irlCurve.__len__()
                k=l=0
                while ((l < lStop) and (k < kStop)):
                    lastK=k
                    lastL=l
                    #TimeStamp Variables
                    A=irlCurve[l][0].__sub__(hgw).__makeInt__()
                    B=selfCurve[k][0].__makeInt__()
                    C=irlCurve[l][0].__add__(hgw).__makeInt__()
                    #FreqStamp Variables
                    X=irlCurve[l][1].__sub__(hfw)
                    Y=selfCurve[k][1]
                    Z=irlCurve[l][1].__add__(hfw)
                    #print "B->",boundPixels,"k,l->",k,l," : ",B," ",A<B,irlCurve[l][0].__makeInt__(),B<C," : ",Y," ",X<Y,irlCurve[l][1],Y<Z
                    if (A < B <= C):
                        k=k+1
                        if (X < Y <= Z):
                            boundPixels=boundPixels+1
                        if (l <lStop-1):
                            if (irlCurve[l][0].__makeInt__() == irlCurve[l+1][0].__makeInt__()):
                                l=l+1
                    elif (B < A):
                        k=k+1
                    else:
                        l=l+1
                fracFound=float(boundPixels).__div__(self.curves[i].length)
                #print "Matching pixels found:",boundPixels," of ",self.curves[i].length," Frac found :",fracFound,"Comparing CurveIDs:",self.curves[i].printCurveID(),iRL.curves[j].printCurveID()
                if (boundPixels == self.curves[i].length):
                    #We remove that candidate from the list
                    candidateToRemove.append(i)
        #Now just copy the candidates which were not found in the clobber
        #listing.
        resultList=copy.deepcopy(self)
        resultList.filename=["ClobberedList"]
        resultList.curves=[]
        for index in range(0,self.curves.__len__()):
            if not candidateToRemove.__contains__(index):
                resultList.curves.append(self.curves[index])
        resultList.totalCount=resultList.curves.__len__()
        #Send back the results
        #print "Keeping :",resultList.totalCount," Removing :",candidateToRemove.__len__(),"Total Candidates Checked :",self.totalCount
        return resultList
    #End clusterClobberWith method

    def candidateStats(self):
        """
        This method calculates stats for the candidate list such as number
        of total curves.  The average curve power and length found and
        their standard deviations.
        Output is a list like
        [startT,stopT,curveCount,AvgP,StdevP,AvgL,StdevL]
        """
        Lsum=float(0)
        Psum=float(0)
        LsumSqr=float(0)
        PsumSqr=float(0)
        curveCount=[]
        gpsStamps=[]
        for entry in self.curves:
            gpsStamps.append(entry.printStartGPS())
            curveCount.append(entry.printCurveID())
            Lsum=Lsum.__add__(float(entry.length))
            LsumSqr=LsumSqr.__add__(float(entry.length).__pow__(2))
            Psum=Psum.__add__(float(entry.power))
            PsumSqr=PsumSqr.__add__(float(entry.power).__pow__(2))
        if (curveCount.__len__() < 1):
            print "Stats can not be performed on candidateless candidateList instance!."
            return []
        gpsStamps.sort()
        meanL=meanP=float(0)
        varL=varP=float(0)
        meanL=Lsum.__div__(curveCount.__len__())
        varL=LsumSqr.__div__(curveCount.__len__())-meanL
        meanP=Psum.__div__(curveCount.__len__())
        varP=PsumSqr.__div__(curveCount.__len__())-meanP
        if gpsStamps.__len__() == 0:
            startT=0
            stopT=0
        else:
            startT=gpsStamps[0]
            stopT=gpsStamps[gpsStamps.__len__()-1]
        statList=[startT,stopT,curveCount.__len__(),meanP,varP,meanL,varL]
        return statList
    #End candidateStats method

    def applyNewThresholds(self,powerEq,pVal,conjunction,lengthEq,lVal):
        """
        Parse these three options to determine a course of action for
        additional thresholds to apply to the candidate lists. Valid options
        are:
        Power Equalities allowed:
        P>float
        P<float
        P=float
        Length Equalities allowed:
        L>integer
        L<integer
        L=integer
        Conjunctions allowed:
        AND  ... Power and Length
        OR   ... Power or Length
        AND! ... A and NOT B
        !AND ... NOT A and B
        myObject.applyNewThresholds('>',3.4,'and','=',8)
        Which means -> ((Power > 3.4) && (Length == 8))
        This method will return an instance of the candidate list which
        fufills the requested criteria.
        NOTE:method __sec2pixel__() is of some use to help turn real units
        seconds to pixels long.
        """
        peq=str(powerEq).lower().lstrip().rstrip()
        pVal=float(pVal)
        Conj=str(conjunction).lower().lstrip().rstrip()
        leq=str(lengthEq).lower().lstrip().rstrip()
        lVal=int(lVal)
        testExpression=peq+Conj+leq
        #>and> ### Power>pVal AND Length>lVal
        newCurveList=[]
        if testExpression == '>and>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower>pVal)and(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<and>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower<pVal)and(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<and<':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower<pVal)and(trackLength<lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '=and=':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower==pVal)and(trackLength==lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '>or>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower>pVal)or(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<or>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower<pVal)or(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<or>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower<pVal)or(trackLength<lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '=or=':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower==pVal)or(trackLength==lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '>and!>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower>pVal)and not(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<and!>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower<pVal)and not(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<and!<':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower<pVal)and not(trackLength<lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '=and!=':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if ((trackPower==pVal)and not(trackLength==lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '>!and>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if (not(trackPower>pVal)and(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<!and>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if (not(trackPower<pVal)and(trackLength>lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '<!and>':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if (not(trackPower<pVal)and(trackLength<lVal)):
                    newCurveList.append(lineInfo)
        elif testExpression == '=!and=':
            for lineInfo in self.curves:
                curveID,trackLength,trackPower=lineInfo.getKurveHeader()
                if (not(trackPower==pVal)and(trackLength==lVal)):
                    newCurveList.append(lineInfo)
        else:
            print 'Error parsing the requested test equalities:',testExpression
            os.abort()
        outputList=self
        outputList.totalCount=newCurveList.__len__()
        outputList.curves=copy.deepcopy(newCurveList)
        return outputList
    #End applyNewThresholds method

#End candidateList class

#Misc methods

def generateFileList(inputTXT):
    """
    This method checks the input string to see if it is
    1) a single file
    2) a list of files
    3) a path
    then we create a list object of files to process
    """
    print "generating file list for :",str(inputTXT)
    absPathFilename=os.path.abspath(inputTXT)
    objList=[]
    dirnameFilename=str(absPathFilename)
    basenameFilename=''
    extensionFilename=''
    if absPathFilename.__contains__("*"):
        dirnameFilename=os.path.dirname(absPathFilename)
        basenameFilename=os.path.basename(absPathFilename)
        extensionFilename=str(str(basenameFilename).split('.')[1])
        print "Found wildcard splitting path into:"
        print dirnameFilename
        print basenameFilename
        print extensionFilename
    if os.path.isdir(dirnameFilename):
        #Create the listing of files to process
        objList=[]
        fileList=[]
        fileList=os.listdir(dirnameFilename)
        for entry in fileList:
            if (str(entry).__contains__(extensionFilename)):
                objList.append(dirnameFilename+'/'+entry)
            #elif (str(entry).__contains__(dirnameFilename)):
            #    objList.append(entry)
    elif os.path.isfile(dirnameFilename):
            #Read this list in from the file
            #Could possibly be just the name of single candidate file!
            objList=[]
            fp=open(dirnameFilename,'r')
            objList=fp.readlines()
            fp.close()
            if str(objList[0]).__contains__('#'):
                #This is a single candidate file specified
                objList=[dirnameFilename]
    else:
        print "Error with getting candidate information from: ",absPathFilename
        objList=[]
    print "file list entries found :",objList.__len__()
    return objList
