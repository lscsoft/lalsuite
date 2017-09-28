#!/usr/bin/env python

# import timeit
# import pandas as pd
from numpy import *
# from pylab import *
import numpy as np
import sys

# # define function
def findCombs2017(threshold,outputName,nround,data):

# findCombs2017(threshold,outputname,nround,data)
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
	elif type(data) is np.ndarray:
		xIn = data
		freq = xIn[:,0]
		power = xIn[:,1]
	else:
		print('Incorrect data type inserted')
		return

	# define threshold value
	sortedPower = np.sort(power)
	if threshold < 1:
		threshInd = np.floor(threshold*len(power))
		thresh = sortedPower[threshInd]
	else:
		thresh = threshold

	epsilon = 0.01				# what is "close enough" to be a comb? ...within epsilon
	f = freq[power >= thresh]	# get only strong frequencies
	f = np.around(f,nround)		# round off to nround decimal points
	print(f)
	f = np.unique(f)			# return only sorted unique values
	nStrong = f.size			# number of strong frequencies


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
			if abs(deltaF_B - deltaF_A) < epsilon: # Are they the same delta f? If so, its a comb!
				if numCombs == 0:
					deltaF = deltaF_A
					comb = np.array([freq0, freq1, freq2])
					numCombs = 1
				else:
					deltaF = np.append(deltaF,deltaF_A)
					comb = np.vstack((comb,[freq0, freq1, freq2]))
					numCombs += 1


	#####################################################
	## Combine three teeth combs with the same delta f ##
	## into one comb with > 3 teeth                    ##
	#####################################################


	deltaF = np.around(deltaF,nround) # Need to round again to account for loss of floating point precision in subtraction
	uniqueDeltaF = np.unique(deltaF) # only need one comb per delta f
	uniqueDeltaF = uniqueDeltaF[np.nonzero(uniqueDeltaF)] # no zero delta fs allowed
	numCombs = len(uniqueDeltaF)

	combList = []
	for combInd in range(numCombs):
		combStorage = comb[abs(deltaF - uniqueDeltaF[combInd]) < epsilon,:]
		numElements = np.size(combStorage)
		combStorage = np.reshape(combStorage,numElements)
		combStorage = np.unique(combStorage)
		combList.append(combStorage)


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
		numTeeth = np.size(combList[combInd])
		# Find highest possible integer multiple of this deltaF that could give a repeat comb
		if numTeeth % 2 == 1: # if odd number of teeth, '%' is modulo
			numRepeatCount = (numTeeth - 1) / 2
		else: # if even number of teeth
			numRepeatCount = (numTeeth - 2) / 2
		# Eliminate all integer multiples that repeat
		for intMult in range(2,numRepeatCount + 1):
			if sum( abs(uniqueDeltaF - intMult*uniqueDeltaF[combInd]) < epsilon ) > 0 : # If there is a repeat comb
				repeatInd = np.where(abs( uniqueDeltaF - intMult*uniqueDeltaF[combInd]) < epsilon) # Find repeat comb index
				repeatInd = int(repeatInd[0]) # make that index an integer
				# combList[repeatInd] = set( combList[repeatInd] ) - set( combList[combInd] ) # set subtract out repeated teeth
				combList[repeatInd] = np.setdiff1d( combList[repeatInd], combList[combInd] ) # set subtract out repeated teeth

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

	fh.write("Combs given as delta f: teeth frequencies ...\n\n")

	for i in range(numCombs):
		if len(combList[i]) != 0:
			outputText = str(uniqueDeltaF[i]) + ": " + str(combList[i])[1:-1] + "\n\n"
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
	findCombs2017(threshold,outputName,nround,data)