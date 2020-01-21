#!/usr/bin/env python

# import timeit
# import pandas as pd
from numpy import *
# from pylab import *
import numpy as np
import sys

# For teeth thicker than epsilon Hz, pick the tooth with the maximum snr.
def uniqueTeeth(epsilon,f,snr):
        indList = [] # List of indices to keep
        thisToothIndices = np.array([0],dtype=np.int) # List of indiced belong to a tooth.
        thisToothSNRs = np.array([snr[0]],dtype=np.int) # List of snrs belong to a tooth.
        for i in range(0, len(f) - 1):
                j = i + 1
                # Use < 2.0*epsilon to avoid round off error to find lines <= epsilon of each other.
                if abs( f[j] - f[i] ) < 2.0*epsilon and j < len(f):
                       thisToothIndices = np.append(thisToothIndices,j)
                       thisToothSNRs = np.append(thisToothSNRs,snr[j])
                else:
                       thisInd = thisToothIndices[np.argmax(thisToothSNRs)]
                       indList = np.append(indList,int(np.floor(thisInd)))
                       # Initialize for next tooth
                       thisToothIndices = np.array([j],dtype=np.int) # List of indiced belong to a tooth.
                       thisToothSNRs = np.array([snr[j]],dtype=np.int) # List of snrs belong to a tooth.
        # We are at the end of the array; add the index of the last tooth.
        thisInd = thisToothIndices[np.argmax(thisToothSNRs)]
        indList = np.append(indList,int(np.floor(thisInd)))
        return indList

# Run this function after the np.unique function.
# This will give a "within epsilon" less strict criteria for what is unique.
# For example, np.unique sees 58.23 Hz and 58.24 Hz as different teeth.
# We want to call them the same tooth.
# Also, after merging combs, can have teeth separated by less than thisDeltaF. Remove those teeth too.
def uniqueWithinEpsilon(epsilon, thisDeltaF, dataList):
    diffList = 1
    for i in range(0, len(dataList) - 1):
        j = i + 1
        #while abs( dataList[j] - dataList[i] ) <= epsilon and j < len(dataList):
        if abs( dataList[j] - dataList[i] ) < 2.0*epsilon or ( thisDeltaF - (dataList[j] - dataList[i]) ) > 2.0*epsilon and j < len(dataList):
            diffList = np.append(diffList,-99999)
            #dataList[j] = -99999
            #j += 1
        else:
            diffList = np.append(diffList,0)
    #dataList = dataList[ dataList != -99999 ]
    #return dataList
    indList = (diffList == 0)
    return indList

# # define function
def findCombs(threshold,outputName,nround,data):

# findCombs(threshold,outputname,nround,data)
# Author: Joe Milliano, STAR Fellow Summer 2017
# Author contact: joe.milliano at gmail dot com
#
# If threshold < 1 then 1-threshold precent of the loudest powers are searched for combs; e.g., threshold=.999 using the top .0001 power values
# If threshold >=1 then po > threshold are serached for combs.
# The outputname get _combs.txt appended to it.
# Give nround equal to, e.g., 2, rounds to frequencies to 2 decimal place.
# If data is a string, it will read the text file with the inputted name and load the data
# If data is a numpy array, it will simply use the data that is there. This method is considerably faster, since loadtxt() is slooooowwwwww


    ########################
    ## Load and sort data ##
    ########################

    np.set_printoptions(formatter={'float': lambda x: "{0:0.2f}".format(x)})

    if type(data) is str:
        xIn = loadtxt(data, skiprows = 1)
        freq = xIn[:,0]
        power = xIn[:,1]
        snr = xIn[:,2]
    elif type(data) is np.ndarray:
        xIn = data
        freq = xIn[:,0]
        power = xIn[:,1]
        snr = xIn[:,2]
    else:
        print('Incorrect data type inserted')
        return

    # define threshold value
    sortedSNR = np.sort(snr)
    if threshold < 1:
        threshInd = int(np.floor(threshold*len(snr)))
        thresh = sortedSNR[threshInd]
    else:
        thresh = threshold

    epsilon = 10 ** (-nround)  # what is "close enough" to be a comb? ...within epsilon
    #epsilon = 0.0             # what is "close enough" to be a comb? ...within epsilon
    f = freq[snr >= thresh]    # get only strong frequencies
    snr = snr[snr >= thresh]   # grab the SNR of the power of the strong frequencies
    f = np.around(f,nround)    # round off to nround decimal points
    #[f, fIdx] = np.unique(f,return_index=True) # Use uniqueTeeth routine below
    #snr = snr[fIdx]
    #print f,len(f)
    #print snr
    # For teeth thicker than epsilon Hz, pick the tooth with the maximum snr.
    uTind = uniqueTeeth(epsilon,f,snr)
    #print uTind
    tmpf = []
    tmpsnr = []
    for i in range(len(uTind)):
        tmpf.append(f[int(uTind[i])])
        tmpsnr.append(snr[int(uTind[i])])
    f = tmpf
    f = np.around(f,nround)     # round off to nround decimal points
    snr = tmpsnr
    snr = np.around(snr,nround) # round off to nround decimal points
    #print f,len(f)
    #print snr
    nStrong = f.size   # number of strong frequencies


    #####################################
    ## Find all combs with three teeth ##
    #####################################


    numCombs = 0
    deltaF = 0
    for ind0 in range(len(f) - 2): # scan all freqs
        freq0 = f[ind0]
        for ind1 in range(ind0 + 1, len(f) - 1 ): # scan all freqs after freq0 (avoids double counting)
            freq1 = f[ind1]
            deltaF_A = freq1 - freq0 # first deltaF
            ind2 = ind1 + 1
            freq2 = f[ind2]
            deltaF_B = freq2 - freq1
            # Get deltaF_B as close to deltaF_A as possible
            if deltaF_A - deltaF_B > epsilon and ind1 + ind2 < nStrong:
                while deltaF_A - deltaF_B > epsilon and ind1 + ind2 <= nStrong:
                    ind2 += 1
                    freq2 = f[ind2]
                    deltaF_B = freq2 - freq1
            if abs(deltaF_B - deltaF_A) <= epsilon: # Are they the same delta f? If so, its a comb!
                if numCombs == 0:
                    deltaF = deltaF_A
                    comb = np.array([freq0, freq1, freq2])
                    indArray = np.array([ind0, ind1, ind2])
                    numCombs = 1
                else:
                    deltaF = np.append(deltaF,deltaF_A)
                    comb = np.vstack((comb,[freq0, freq1, freq2]))
                    indArray = np.vstack((indArray,[ind0, ind1, ind2]))
                    numCombs += 1


    #####################################################
    ## Combine three teeth combs with the same delta f ##
    ## into one comb with > 3 teeth                    ##
    #####################################################


    deltaF = np.around(deltaF,nround) # Need to round again to account for loss of floating point precision in subtraction
    uniqueDeltaF = np.unique(deltaF) # only need one comb per delta f
    # snr = snr[uniqueFIdx]
    # snr = snr[np.nonzero(uniqueDeltaF)]
    uniqueDeltaF = uniqueDeltaF[np.nonzero(uniqueDeltaF)] # no zero delta fs allowed
    numCombs = len(uniqueDeltaF)

    combList = []
    snrList = []
    for combInd in range(numCombs):
        rowInd = abs( deltaF - uniqueDeltaF[combInd] ) <= epsilon # Which rows hold teeth of the same comb?
        if numCombs == 1:
            combStorage = comb
            indStorage = indArray
        else:
            combStorage = comb[rowInd,:]        # Grab those teeth
            indStorage = indArray[rowInd,:]     # Also grab their indices from the original freq vector...this is used to get SNR of the tooth signal
            numElements = np.size(combStorage)  # How many teeth are in this comb?
            combStorage = np.reshape(combStorage,numElements) # Turn teeth matrix into a vector
            indStorage = np.reshape(indStorage, numElements)  # Turn index matrix into a vector
            [combStorage, toothIdx] = np.unique(combStorage,return_index=True) # Eliminate repeat teeth
            indStorage = indStorage[toothIdx]   # Eliminate indices of repeat teeth
            #combStorage = [combStorage[tInd] for tInd in sorted(toothIdx)] # Put the unique teeth in the order found
            #combStorage = np.around(combStorage,nround) # Need to round again.
        snrStorage = snr[indStorage]    # get SNR of power signal at that tooth
        combList.append(combStorage)
        snrList.append(snrStorage)

    for i in range (len(combList)):
        #combList[i] = uniqueWithinEpsilon(epsilon, combList[i])
        thisComb = combList[i]
        thisDeltaF = uniqueDeltaF[i]
        thisIndList = uniqueWithinEpsilon(epsilon,thisDeltaF,thisComb)
        combList[i] = thisComb[thisIndList]

    ########################################
    ## Eliminate higher freq repeat teeth ##
    ########################################


    ## KEY INSIGHT:
    ## Let's say a particular deltaF has N teeth in a comb.
    ## Then an integer multiple of deltaF can repeat some of these N combs.
    ## The largest integer n that can have repeat combs is
    ## the largest n such that 2n < N-1.

    ## COUNTING REPEAT COMB STEPS:
    ## 1. Go to integer multiple of deltaF, say intMult*deltaF
    ## 2. Find and eliminate all matching teeth for that frequency (the final set subtraction)
    ## 3. Repeat for all possible integer multiples of deltaF (inner loop)
    ## 4. Repeat for all deltaF (outer loop)

    for combInd in range(numCombs):
        thisDeltaF = uniqueDeltaF[combInd]
        # Loop through the combs with smaller deltaF, and check if this deltaF is a harmonic of the smaller deltaF.
        for smallerDeltaFInd in range(combInd):
            smallerDeltaF = uniqueDeltaF[smallerDeltaFInd]
            if abs(thisDeltaF - smallerDeltaF*around((thisDeltaF/smallerDeltaF),0)) <= 2.0*epsilon:
                combList[combInd] = np.setdiff1d( combList[combInd], combList[smallerDeltaFInd] ) # subtract out repeated teeth in harmonic of a comb

    ###################################################
    ## If comb frequencies are within epsilon Hz,       ##
    ## and if they have the same teeth,              ##
    ## then call them the same comb (rounding error) ##
    ###################################################

    deleteInd = [] # Index of rows we will have to delete
    for i in range(len(uniqueDeltaF)-1) :
        if len(combList[i]) > len(combList[i+1]) : # Keep as many teeth as possible
            ind0 = i + 1
            ind1 = i
        else:
            ind0 = i
            ind1 = i + 1
        if abs( uniqueDeltaF[ind0] - uniqueDeltaF[ind1] ) < 2*epsilon :
            combList[ind0] = np.setdiff1d( combList[ind0], combList[ind1] )
            combList[ind0] = np.setdiff1d( combList[ind0] - epsilon, combList[ind1] )
            combList[ind0] = np.setdiff1d( combList[ind0] + epsilon, combList[ind1] )
            if len( combList[ind0] ) < 3 :
                deleteInd = np.append(deleteInd,ind0)

    if len(deleteInd) > 0 :
        deleteInd = deleteInd.astype(int) # Make indices into integers
        uniqueDeltaF = np.delete(uniqueDeltaF,deleteInd)
        combList = np.delete(combList,deleteInd)
        numCombs = len(uniqueDeltaF)




    #######################
    ## Write Output File ##
    #######################


    fh = open(outputName + "_combs.txt","w")

    fh.write("Combs given as delta f (median SNR): teeth frequencies ...\n\n")

    for i in range(numCombs):
        if len(combList[i]) > 2:
            medianSNR = np.median(snrList[i])
            medianSNR = np.around(medianSNR,nround)
            # If there is a comb within a comb with the same spacing but shifted back in frequency, add a comma to separate them in the list.
            thisComb = combList[i]
            thisCombStr = str(thisComb[0])
            for j in range(1,len(thisComb)):
                if thisComb[j] < thisComb[j-1]:
                    thisCombStr = thisCombStr + ', ' +  str(thisComb[j])
                else:
                    thisCombStr = thisCombStr + ' ' +  str(thisComb[j])
            #outputText = str(uniqueDeltaF[i]) + " Hz (" + str(medianSNR) + "): " + str(combList[i])[1:-1] + "\n\n"
            outputText = str(uniqueDeltaF[i]) + " Hz (" + str(medianSNR) + "): " + thisCombStr + "\n\n"
            fh.write(outputText)

    fh.close()


######################################################
## If this is the main function being run, then     ##
## get arguments from command line and run function ##
######################################################


if __name__ == '__main__':
    # get arguments
    threshold = sys.argv[1]
    outputName = sys.argv[2]
    nround = sys.argv[3]
    data = sys.argv[4]

    # turn relevant strings into floats, ints
    if isinstance(threshold,str):
        threshold = float(sys.argv[1])
    if isinstance(nround,str):
        nround = int(sys.argv[3])

    # run function
    findCombs(threshold,outputName,nround,data)
