#!/usr/bin/env python

from numpy import *
import matplotlib
matplotlib.use('Agg')
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from pylab import *
import scipy as scipy
import sys
import os
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
import math as math
import numpy as np


def plotSpecAvgOutput(filename,outputFileName,chanName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coinDF,referenceFile):

# fileName       -- the name of the file with data to input; this file is the output from spec_avg.
# outputFileName -- base name of the output spectrogram files; will append .pdf and .png to this name.
# chanName       -- the name of the channel used to generate the spectrogram data.
# effTBase       -- the effective time base line: 1/effTBase gives the frequency resolution of the plot.
#                   so that a effTBase of 10 seconds means deltaF = 0.1 Hz.
# deltaFTicks    -- the change in frequency between major tick marks (e.g., 5 Hz).
# taveFlag       -- if > 0 then produce StackSlide style time average output without sliding.
#                   The spectrograms is given as subplot(2,1,1) and this spectrum as subplot(2,1,2).
#                   Also, .txt is appended to the plot name, and f versus the normalized averaged
#                   power is output into this text file.
# effTBaseFull   -- 1/effTBaseFull gives the frequency resolution of the average normalized spectra plots. This value is the timbaseline of the sfts.
# thresholdSNR   -- if > 0 then look for coincident lines with the referenceFile spectra above this threshold.
# coinDF         -- window in frequency to use when looking for coincident lines.
# referenceFile  -- base name of the reference file output by spec_avg; will append _timeaverage to this name.

# Convert relevant strings to numbers.
    if isinstance(effTBase,str):
        effTBase = float(effTBase)
    if isinstance(deltaFTicks,str):
        deltaFTicks = int(deltaFTicks)
    if isinstance(taveFlag,str):
        taveFlag = float(taveFlag)
    if isinstance(effTBaseFull,str):
        effTBaseFull = float(effTBaseFull)
    if isinstance(thresholdSNR,str):
        thresholdSNR = float(thresholdSNR)
    if isinstance(coinDF,str):
        coinDF = float(coinDF)


    xIn=loadtxt(filename)        #Loads contents of the file into xIn
    xlen=len(xIn[:,0])
    ylen=len(xIn[0,:])

    y=xIn.transpose()[ ::-1,:]   #Transpose and flip rows upside-down
    lst=filename.split("_")      #Breaks up filename into list of elements seperated by "_"

    tEnd = lst.pop()        # end time
    tStart = lst.pop()      # start time
    ifo = lst.pop()         # ifo
    fEnd = lst.pop()        # end frequency
    fStart = lst.pop()      # start frequency 

    y_temp1 = y             #Create an array of the spectrogram data without the segments taken out when a segment file is used
    yzeros=list(y.sum(axis=0))
    for i in reversed(range(len(yzeros))):
        if yzeros[i] == 0:
            y_temp1 = scipy.delete(y_temp1,i,1)

    y_temp2 = []
    xlen2=len(y_temp1[:,0])
    ylen2=len(y_temp1[0,:])
    for i in range(xlen2):          #puts y_temp1 in a list column by column
        for j in range(ylen2):
            y_temp2.append(y_temp1[i,j])


    filename_numbins = filename + '_numbins'
    NumBinsAvg=loadtxt(filename_numbins)

    mediany = median(y_temp2)
    cutoffval = mediany + (5*(mediany/sqrt(NumBinsAvg)))

#Replace all of the values greater that cut-off value with cutoffval
    for i in range(len(y[:,0])):
        for j in range(len(y[0,:])):
            if y[i,j] >= cutoffval:
                y[i,j] = cutoffval

# Get dates, etc.
    filename_date = filename + '_date'
    dates = loadtxt(filename_date)

#Create strings of start and ends times in 'YYYY/MM/DD HH:MM:SS' format
    startTime = str(int(dates[0,1]))+'/'+str(int(dates[0,2]))+'/'+str(int(dates[0,3]))+' '+str(int(dates[0,4]))+':'+str(int(dates[0,5]))+':'+str(int(dates[0,6]))
    endTime = str(int(dates[-1,1]))+'/'+str(int(dates[-1,2]))+'/'+str(int(dates[-1,3]))+' '+str(int(dates[-1,4]))+':'+str(int(dates[-1,5]))+':'+str(int(dates[-1,6]))

#---------------------------------
#Plot spectrogram
#---------------------------------

    figure(1)

#Modify jet colormap
    ymax = y.max()
    ymin = y_temp1.min()                     #min can be zero
    yrange = ymax - ymin
    cymin = ymin - (yrange/254)              #Create a range of 254 steps instead of 256
    jetcmap = get_cmap('jet')
    cvalues = jetcmap(arange(256))
    cvalues[0] = [0.85,0.85,0.85,1.]         #Define colormap so that values under minimum will show up as grey
    mod_cmap = cm.colors.LinearSegmentedColormap.from_list("newjet", cvalues)
#Plot spectrogram
    imshow(y,aspect='auto',cmap=mod_cmap,interpolation='nearest',vmin=cymin,vmax=ymax)

#    jetcmap.set_under('0.85')
#    imshow(y,aspect='auto',cmap=jetcmap)
    cbar=colorbar(format='%d')             #Show colorbar on spectrogram
#    clim(ymin,ymax)
    cbar.set_label('Channel units/root Hz')  #Colorbar label
    filename_png = filename + '.png'
    
#y-axis
    ylabel('Frequency [Hz]',fontsize=10)
    fRange = int(float(fEnd)) - int(float(fStart))
    VecFLabels = range(int(float(fStart)),int(float(fEnd))+1,deltaFTicks)
    VecFLabels.reverse()
    yticks(range(0,int((fRange + 1)*effTBase),int(deltaFTicks*effTBase)),VecFLabels,fontsize=9)

#x-axis
    xlabel('Time in days since start date',fontsize=10)
#    deltaxticks = int(round(xlen/10))
    deltaxticks = 10
    if deltaxticks < 1:
        deltaxticks = 1
    xtcks = linspace(-.5,xlen-.5,deltaxticks)      #matplotlib on cluster seems to start counting tick locations at -.5 insteasd of zero
    xtcks = xtcks.tolist()
    xnumticks=len(xtcks)
    tRange = float(tEnd)-float(tStart)
    xlabs=linspace(0,tRange,xnumticks)
    xlabs = xlabs.tolist()
    for i in range(len(xlabs)):
        day=xlabs[i]/86400
        xlabs[i] = "%.2f" %day
    xticks(xtcks,xlabs,fontsize = 9)

#Title

    titleString = 'Spectrogram for %s;  %s  to  %s UTC.' %(chanName,startTime,endTime)
    title(titleString,fontsize = 9)

#Save spectrogram
    savefig(filename_png)
    close()

#---------------------------------
#Plot normalized average power
#---------------------------------

    figure(2)

    timeaverageFileName = filename + '_timeaverage'
    data = loadtxt(timeaverageFileName)       #Load data for normalized average power plot
    fk = data[:,0]                            #Frequency data
    fklist = fk.tolist()
    xout = data[:,1]                          #Normalized average power data
    xoutlist = xout.tolist()
    snr = []
    for i in range(len(xout)):                #Calculate SNR using number of sft's used (don't include unused timesteps from seg. list)
        sr = (xout[i]-1)*math.sqrt(ylen2)
        snr.append(sr)
    data2=zeros((len(fk),3))
    for i in range(len(xout)):
        data2[i,0] = data[i,0]
        data2[i,1] = data[i,1]
        data2[i,2] = snr[i]

    outputTextFile = filename + '.txt'
    outputSortedTextFile = outputFileName + '_sorted.txt'

    ind=lexsort((data2[:,0],data2[:,1],data2[:,2]))         #Get indexes for sorted data
    data3=data2[ind]                            #Make sorted 2d array
    data4=data3[::-1]                          #Flip array
    fSorted = data4[:,0]
    xoutSorted = data4[:,1]
    snrSorted = data4[:,2]

    f_header = open(outputTextFile,'a')
    f_header.write('Freq.     Normalized Avg Power   SNR \n')
    savetxt(f_header,transpose(array((fk,xout,snr))),fmt="%.6f %.4f %.4f")
    f_headers = open(outputSortedTextFile,'a')
    f_headers.write('Freq.     Normalized Avg Power   SNR\n')
    savetxt(f_headers,transpose(array((fSorted,xoutSorted,snrSorted))),fmt="%.6f %.4f %.4f")
    kmax = len(xout)
    stdev_xout = std(xout)
    meanval_xout = mean(xout)

# Read in timestamps file to find the number of SFTs used:
    timestampFileName = filename + '_timestamps'
    tmp = loadtxt(timestampFileName)
    ntmp = tmp[:,0]
    ttmp = tmp[:,1]

# Computed expected 5 sigma cutoff for gaussian noise:
    cutoffmax = 1.0 + 5.0/sqrt(len(ntmp))
# avoid too small of a maximum cutoff:
    if cutoffmax < 4:
        cutoffmax = 4

# Compute cutoff from data, but do not exceed cutoffmax:
    cutoff = meanval_xout+(5*stdev_xout)
    if cutoff > cutoffmax:
        cutoff = cutoffmax

    for i in range(len(xoutlist)):
        if xoutlist[i] >= cutoff:
            xoutlist[i] = cutoff

#Plot normalized average power vs. frequency
    plot(fklist,xoutlist,'-',linewidth=0.6)
#    axis([float(fStart), float(fEnd), 0., 4.])
#Take care of labels and axes
    titleString = 'Spectrum for %s; averaged over %s to %s UTC.' %(chanName,startTime,endTime)
    title(titleString,fontsize = 9)
    ylabel('Normalized Average Power',fontsize = 10)
    yticks(fontsize = 9)
#    majorLocator   = MultipleLocator(100)
#    majorFormatter = FormatStrFormatter('%d')
#    minorLocator   = MultipleLocator(50)
#    xaxis.set_major_locator(majorLocator)
#    xaxis.set_major_formatter(majorFormatter)
#    xaxis.set_minor_locator(minorLocator)
    plt.locator_params(axis = 'x', nbins = 11)
    xlabel('Frequency (Hz)',fontsize = 10)
    xticks(fontsize = 9)
#    XFLabels = range(int(float(fStart)),int(float(fEnd))+1,50)
#    xticks(range(0,(fRange + 1),50),XFLabels)
#    figure(figsize=(20, 20))
    filename_png2 = filename + '_2.png'

#Save normalized average power plot
    savefig(filename_png2)
    close()
#    savefig(filename_png2, dpi=900)
    if thresholdSNR > 0:
        timeaverageFileNameRef = referenceFile + '_timeaverage'
        refData = loadtxt(timeaverageFileNameRef)       #Load reference data
        fRef = refData[:,0]
        xRef = refData[:,1]

        outputTextFileLines = outputFileName + '_coincident_lines.txt'
        fid3 = open(outputTextFileLines,'w')
        fid3.write('\n       COINCIDENT LINES       \n')
        fid3.write('\n')
        fid3.write('             INPUT                             REFERENCE          \n')
        fid3.write(' Freq. (Hz)   Power      SNR       Freq. (Hz)   Power      SNR    \n')
        fid3.write('\n')

        outputTextFileNewLines = outputFileName + '_new_lines.txt'
        fid4 = open(outputTextFileNewLines,'w')
        fid4.write('\n       NEW LINES       \n')
        fid4.write('\n')
        fid4.write('             INPUT                 \n')
        fid4.write(' Freq. (Hz)   Power      SNR       \n')
        fid4.write('\n')

        outputTextFileOldLines = outputFileName + '_old_lines.txt'
        fid5 = open(outputTextFileOldLines,'w')
        fid5.write('\n       OLD LINES       \n')
        fid5.write('\n')
        fid5.write('             INPUT                 \n')
        fid5.write(' Freq. (Hz)   Power      SNR       \n')
        fid5.write('\n')

        xoutMean = mean(xout)
        oneOverxoutSTD = 1.0/std(xout)
        xRefMean = mean(xRef)
        oneOverxRefSTD = 1.0/std(xRef)
        SNRout = (xout - xoutMean)*oneOverxoutSTD
        SNRRef = (xRef - xRefMean)*oneOverxRefSTD
        lengthSNRout = SNRout.size
        skip = 0
        iMax = 0
        coincidenceBins = math.ceil(coinDF*effTBaseFull)
        for j in range(lengthSNRout):
            if skip == 0:
                jMin = j - coincidenceBins
                if jMin < 1:
                    jMin = 1
                jMax = j + coincidenceBins
                if jMax > lengthSNRout:
                    jMax = lengthSNRout
                shortRange = range(int(jMin),int(jMax+1))
                indexRangeList1 = []
                for i in shortRange:
                    indexRangeList1.append(SNRout[i])
                SNRoutmax = 0
                iMaxout = 0
                if len(indexRangeList1) > 0:
                    SNRoutmax = max(indexRangeList1)
                    iMaxout = indexRangeList1.index(max(indexRangeList1))
                    iMaxout = jMin + iMaxout - 1
                else:
                    skip = 1
                indexRangeList2 = []
                for i in shortRange:
                    indexRangeList2.append(SNRRef[i])
                SNRRefmax = 0
                iMaxRef = 0
                if len(indexRangeList2) > 0:
                    SNRRefmax = max(indexRangeList2)
                    iMaxRef = indexRangeList2.index(max(indexRangeList2))
                    iMaxRef = jMin + iMaxRef - 1
                else:
                    skip = 1
                if SNRoutmax >= thresholdSNR and SNRRefmax >= thresholdSNR:
                    skip = 1
                    fid3.write(' %11.6f  %9.4f  %7.2f    %11.6f  %9.4f  %7.2f\n' % (fk[iMaxout],xout[iMaxout],SNRoutmax,fRef[iMaxRef],xRef[iMaxRef],SNRRefmax))
                elif SNRoutmax >= thresholdSNR and SNRRefmax < thresholdSNR:
                    skip = 1
                    fid4.write(' %11.6f  %9.4f  %7.2f\n' % (fk[iMaxout],xout[iMaxout],SNRoutmax))
                elif SNRoutmax < thresholdSNR and SNRRefmax >= thresholdSNR:
                    skip = 1
                    fid5.write(' %11.6f  %9.4f  %7.2f\n' % (fk[iMaxout],xout[iMaxout],SNRoutmax))
            else:
                if (j - iMaxRef) > coincidenceBins and (j - iMaxout) > coincidenceBins:
                    skip = 0
        fid3.close()
        fid4.close()
        fid5.close()

from findCombs2017 import findCombs2017

#Load data from command line
filename = sys.argv[1]
outputFileName = sys.argv[2]
chanName = sys.argv[3]
effTBase = sys.argv[4]
deltaFTicks = sys.argv[5]
taveFlag = sys.argv[6]
effTBaseFull = sys.argv[7]
thresholdSNR = sys.argv[8]
coinDF = sys.argv[9]
referenceFile = sys.argv[10]

#Call plotting function and comb finding function
plotSpecAvgOutput(filename,outputFileName,chanName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coinDF,referenceFile)
findCombs2017(0.999,outputFileName,2,filename + "_timeaverage")