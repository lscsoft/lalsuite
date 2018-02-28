#!/usr/bin/env python
# -*- coding: utf-8 -*-
#       cbcBayesCombinePosteriors.py
#
#       Copyright 2016
#       Christopher Berry <christopher.berry@ligo.org>
#       Sebastian Gaebel <sebastian.gaebel@ligo.org>
#
#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 2 of the License, or
#       (at your option) any later version.
#
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#
#       You should have received a copy of the GNU General Public License
#       along with this program; if not, write to the Free Software
#       Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#       MA 02110-1301, USA.

"""
This is used for combining posterior samples together.
"""

#Import standard things
import argparse
try:
    from lalinference import LALInferenceHDF5PosteriorSamplesDatasetName as posterior_grp_name
except:
    posterior_grp_name = "posterior_samples"
from collections import defaultdict
import h5py
import numpy as np

#Set-up commands
parser = argparse.ArgumentParser(description="Combine some posterior samples.")
parser.add_argument("-o", "--output", dest="outfilename", default="combined_posterior",
                  help="Write combined posterior to OUTFILE", metavar="OUTFILE")
parser.add_argument("-p", "--pos", action="append", dest="infilename",
                  help="Combine posteriors from INFILE", metavar="INFILE", required=True)
parser.add_argument("-a", "--all", dest="all",
                  action="store_true", default=False,
                  help="Use all posterior samples (do not weight)")
parser.add_argument("-w", "--weight", action="append", default=[], dest="weightings",
                  help="Weighting for posteriors", type=float)
shuffleGroup = parser.add_mutually_exclusive_group()
shuffleGroup.add_argument("-s", "--shuffle",
                  action="store_true", dest="shuffle", default=True,
                  help="Randomise posterior samples before combination [Default]")
shuffleGroup.add_argument("-ns", "--no-shuffle",
                  action="store_false", dest="shuffle",
                  help="Do not randomise posterior samples before combination")
parser.add_argument("-m", "--mix", dest="mix",
                  action="store_true", default=False,
                  help="Randomise combined samples")
fileGroup = parser.add_mutually_exclusive_group()
fileGroup.add_argument("-t", "--text",
                  action="store_false", dest="hdf", default=False,
                  help="Use ASCII posterior (.dat) files [Default]")
fileGroup.add_argument("-f", "--hdf",
                  action="store_true", dest="hdf",
                  help="Use HDF5 posterior files")


args = parser.parse_args()

#Count arguments
nPos = np.size(args.infilename)
nWeight = np.size(args.weightings)

print args.weightings

#Check sensible combination of arguments
if (nWeight != 0):
    if args.all:
        print "You cannot use all posterior samples and weight them!"
        exit(1)

    if (nWeight != nPos):
        print "Please either specify a weight for each posterior file or none"
        exit(1)
else:
    args.weightings = [1.0] * nPos

#Specify combination ID
combineID = "combined"
if args.all:
    combineID = combineID+"_all"
else:
    combineID = combineID+"_weight_"+'_'.join(map(str, args.weightings))

if args.shuffle:
    combineID = combineID+"_shuffle"
else:
    combineID = combineID+"_noshuffle"

if args.mix:
    combineID = combineID+"_mixed"

print "Combined ID:", combineID

#Initiate lists to hold data 
samples = []
paramsList = []
sizeList = []
metadata = {"lalinference": defaultdict(lambda: [None]*nPos),
            "lalinference/"+combineID: defaultdict(lambda: [None]*nPos),
            "lalinference/"+combineID+"/"+posterior_grp_name: defaultdict(lambda: [None]*nPos)}

#Read in data
for posIndex in range(nPos):
    if args.hdf:
        #HDF5 files with metadata
        with h5py.File(args.infilename[posIndex], "r") as inFile:

            group = inFile["lalinference"]
            for key in group.attrs:
                metadata["lalinference"][key][posIndex] = group.attrs[key]

            run_id = list(group.keys())[0]

            group = group[run_id]

            for key in group.attrs:
                metadata["lalinference/"+combineID][key][posIndex] = group.attrs[key]

            if "combined_run_ids" not in group.attrs:
                metadata["lalinference/"+combineID]["combined_run_ids"][posIndex] = run_id
            elif "recombined_run_ids" not in group.attrs:
                metadata["lalinference/"+combineID]["recombined_run_ids"][posIndex] = run_id
            elif "rerecombined_run_ids" not in group.attrs:
                metadata["lalinference/"+combineID]["rerecombined_run_ids"][posIndex] = run_id
            elif "rererecombined_run_ids" not in group.attrs:
                metadata["lalinference/"+combineID]["rererecombined_run_ids"][posIndex] = True
            else:
                print "Too many combinations to count!"


            group = group[posterior_grp_name]
            for key in group.attrs:
                metadata["lalinference/"+combineID+posterior_grp_name][key][posIndex] = group.attrs[key]

            posDtype = []
            for key in group:
                posDtype.append((key, group[key].dtype))
                shape = group[key].shape

            posData = np.empty(shape, dtype=posDtype)

            for key in group:
                posData[key] = group[key][:]

    else:
        #Standard text file
        posData = np.genfromtxt(args.infilename[posIndex], names=True)

    if (args.shuffle):
        np.random.shuffle(posData)

    samples.append(posData)
    paramsList.append(set(posData.dtype.names))
    sizeList.append(np.size(posData))


#Create intersection
paramsOut = list(set.intersection(*paramsList))

datatypes = samples[0][paramsOut].dtype

#Combine posteriors
if (args.all):
    #Use all samples
    sizeOut = sum(sizeList)
    samplesOut = np.empty(sizeOut, dtype=datatypes)

    indexSize = sizeList
    
else:
    #Weight different posteriors
    fracWeight = np.asarray(args.weightings) / float(sum(args.weightings))

    testNum = fracWeight * float(sum(sizeList))
    minIndex = np.argmin(np.asarray(sizeList) / np.asarray(testNum))
    
    testSize = sizeList[minIndex] / fracWeight[minIndex]

    weightNum = np.around(fracWeight * testSize)
    sizeOut = sum(weightNum)
    samplesOut = np.empty(sizeOut, dtype=datatypes)

    indexSize = weightNum


print "Using number of samples ", indexSize


startIndex = 0
for posIndex in range(0,nPos):
    stopIndex = startIndex + indexSize[posIndex]

    for paramIndex, paramItem in enumerate(paramsOut):
        samplesOut[paramItem][startIndex:stopIndex] = samples[posIndex][paramItem][0:indexSize[posIndex]]  

    startIndex = stopIndex


#Mix samples
if args.mix:
    np.random.shuffle(samplesOut)


#Save output
if args.hdf:
    #HDF5 file with metadata
    with h5py.File(path, "w") as outFile:
        group = outFile.create_group("lalinference")
        group = group.create_group(combineID)
        group = group.create_group(posterior_grp_name)
        for key in samplesOut.dtype.names:
            group.create_dataset(key, data=samplesOut[key], shuffle=True, compression="gzip")

        for level in metadata:
            for key in metadata[level]:
                outFile[level].attrs[key] = metadata[level][key]

else:
    #Standard textt output
    paramHeader = "\t".join(paramsOut)
    np.savetxt(args.outfilename, samplesOut.T, delimiter="\t", header=paramHeader, comments="")


#Done!

