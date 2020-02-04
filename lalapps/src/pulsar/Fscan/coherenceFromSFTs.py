#!/usr/bin/env python

import os
import struct
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')

from matplotlib.pyplot import *

import matplotlib.pyplot as plt

matplotlib.rcParams['text.usetex']=False

#Written by Kara Merfeld, 2018 (In my defense, I was a little first year grad student when I did this :) )

# parseSFT function from John Whelan
def parseSFT(SFTinput):
    """
    This function parses an SFT (short fourier transform) formatted
    according to the SFT specification version 2,
    https://dcc.ligo.org/cgi-bin/private/DocDB/ShowDocument?docid=T040164
    It takes as input either an open (for binary reading) file object
    or a string specifying the path to the SFT file.  It returns the
    SFT data and metadata in a dictionary.  If more than one SFT is in
    a file, it returns a tuple of dictionaries.
    """
    if type(SFTinput) == str:
        SFTfile = open(SFTinput,'rb')
    elif type(SFTinput) == file:
        SFTfile = SFTinput
    else:
        raise TypeError('Argument must be a string (filename) or open file object.')

    sfts = []
    while True: # Loop over SFTs in the file

        header_data = SFTfile.read(48)
        if len(header_data) == 0: # We've reached the end
            break

        # Parse header (48 bytes)
        ( version,
          gps_sec,
          gps_nsec,
          tbase,
          first_frequency_index,
          nsamples,
          crc64,
          detector,
          padding,
          comment_length
          ) = struct.unpack('diidiiQ2s2si', header_data)

        # Check version (and endianness)
        if version != 2.0:
            raise ValueError('Can only parse SFTs of version 2, not %f'
                               % version)

        # TODO: check crc64 checksum

        # Read comment
        comment = struct.unpack( ( '%ds' % comment_length ),
                                 SFTfile.read(comment_length))[0]

        # Read data
        data_stream = struct.unpack( ( '%df' % nsamples*2 ),
                                     SFTfile.read(nsamples*8))

        # Pack real values into complex array
        data = np.array(data_stream[0::2]) + 1.0j * np.array(data_stream[1::2])

        sft = {
            'version' : version,
            'gps_sec': gps_sec,
            'gps_nsec' : gps_nsec,
            'tbase' : tbase,
            'first_frequency_index' : first_frequency_index,
            'nsamples' : nsamples,
            'crc64' : crc64,
            'detector' : detector,
            'padding' : padding,
            'comment_length' : comment_length,
            'comment' : comment,
            'data' : data
            }

        sfts.append(sft)

    if len(sfts) == 1:
        return sfts[0]
    else:
        return tuple(sfts)


def coherenceFromSFTs( pathToSFTsChannA, pathToSFTsChannB, subBand=100):  #The function that generates the coherence

    done = False
    for i in range(11 , len(pathToSFTsChannA)-1):
        if pathToSFTsChannA[-i]=='/' and done == False:
            CA = pathToSFTsChannA[-(i-1):-10]
            done = True

    done = False
    for i in range(11 , len(pathToSFTsChannB)-1):
        if pathToSFTsChannB[-i]=='/' and done == False:
            CB = pathToSFTsChannB[-(i-1):-10]
            done = True

    print('Computing the coherence between:')
    print(CA)
    print(CB)


    #name = pathToSFTsChannA[-20:-10]
    #print(name)

    ListA = sorted(os.listdir(pathToSFTsChannA))
    ListB = sorted(os.listdir(pathToSFTsChannB))

    channelA = []
    channelB = []

    #We create lists of the starting times:
    StartTimesA = []
    StartTimesB = []  #Ultimately we will need to account for channels with different starting frequencies too.


    ##We find the starting frequency

    FFI_A =  parseSFT(pathToSFTsChannA+ListA[0])['first_frequency_index'] #the first frequency indices
    FFI_B =  parseSFT(pathToSFTsChannB+ListB[0])['first_frequency_index']

    TbaseA = parseSFT(pathToSFTsChannA+ListA[0])['tbase']
    TbaseB = parseSFT(pathToSFTsChannB+ListB[0])['tbase']

    #StartFreqA = parseSFT(pathToSFTsChannA+ListA[0])['data'][FFI_A]
    #StartFreqB = parseSFT(pathToSFTsChannB+ListB[0])['data'][FFI_B]

    FminA = FFI_A/TbaseA
    FminB = FFI_B/TbaseB

    Fmin = 0
    if FminA <= FminB:
        Fmin = FminB
    else:
        Fmin = FminA

    ##We find the ending frequency:

    FmaxA = FminA + (parseSFT(pathToSFTsChannA+ListA[0])['nsamples']-1)/TbaseA
    FmaxB = FminB + (parseSFT(pathToSFTsChannB+ListB[0])['nsamples']-1)/TbaseB

    Fmax = 0

    if FmaxA <= FmaxB:
        Fmax = FmaxA
    else:
        Fmax = FmaxB

    #We calculate the beginning and ending bins
    KminA = int((Fmin-FminA)*TbaseA)
    KminB = int((Fmin-FminB)*TbaseB)

    KmaxA = int((Fmax-FminA)*TbaseA)
    KmaxB = int((Fmax-FminB)*TbaseB)

    #print('K values')
    #print(KminA)
    #print(KmaxA)
    #print(KminB)
    #print(KmaxB)

    #A = A[KminA:KmaxA+1]
    #B = B[KminB:KmaxB+1]


    for sft in ListB:
        StartTimesB.append(float(sft[-19:-10]))

    A = 0 * parseSFT(pathToSFTsChannA+ListA[0])['data'] #gotta come back and make sure the lengths work out if they are not equal
    B = 0 * parseSFT(pathToSFTsChannB+ListB[0])['data']
    #numerator = 0 * parseSFT(pathToSFTsChannA+ListA[0])['data']

    #print('lengths of A and B before slice:')
    #print(len(A))
    #print(len(B))


    A = A[KminA:KmaxA+1]
    B = B[KminB:KmaxB+1]

    numerator = A

    #print('lengths')
    #print(len(A))
    #print(len(B))
    #print(len(numerator))


    for sftA in ListA:
        StartTimesA.append(float(sftA[-19:-10])) #This list might not even be necesarry

    nAve = 0
    #Let A be the channel that the thing is getting compared to.
    for Aind in range(0,len(StartTimesA)):
        for Bind in range(0,len(StartTimesB)):
            if StartTimesA[Aind] == StartTimesB[Bind]:

                channelA = parseSFT(pathToSFTsChannA+ListA[Aind])['data'][KminA:KmaxA+1]
                channelB = parseSFT(pathToSFTsChannB+ListB[Bind])['data'][KminB:KmaxB+1]

                A = A + channelA * np.conj(channelA)
                B = B + channelB * np.conj(channelB)
                numerator = numerator + channelA * np.conj(channelB)
                nAve += 1 # Keep track of the number of averages


    numerator = numerator * np.conj(numerator)

    coh = numerator/(A*B)
    coh = np.real_if_close(coh, tol = 10)

    nAve = str(nAve)

    print('Coherence Completed; nAve = %s' % nAve)
    print('Generating plots and files')
    ###

    T = TbaseA
    B = Fmax-Fmin  #I am not understanding what B is... It is the bandwidth of frequencies -- it is
    N = T * B +1 # Number of frequency bins? I think so, because N=len(coh)

    Freq = np.linspace(Fmin,Fmax, N)

    #print('lenth of frequency array and coherence array')
    #print(len(Freq))
    #print(len(coh))

    #Freq = np.linspace(0, N, len(coh))
    #print('value of Freq[180000]')
    #print(Freq[180000])   #this is telling me that at each index is 1 Hz.  No good.

    #Sampling frequency F = 1/1800. So then to get to 100Hz, we'd go through 180000 indices


    #Now we create a loop for the Figures:
    #We want this to run from Fmin to Fmax

    detector = parseSFT(pathToSFTsChannA+ListA[0])['detector']
    detector = str(detector)

    startTime = parseSFT(pathToSFTsChannA+ListA[0])['gps_sec']
    #endTime = startTime + TbaseA *len(ListA)

    startTime = str(startTime)
    #endTime = str(endTime)

    #We create an array for the coherence output:

    Coh_Output = np.concatenate((Freq.reshape(-1,1),coh.reshape(-1,1)), axis = 1)
    #print(Coh_Output.shape)

    #print(Coh_Output[0:3])

    subBand = int(subBand) # Output plots and files for each subBand.

    # All frequencies below, minFreq, maxFreq, and subBand, are coverted to integer indices in the Freq array.
    i = int(0)
    subBand = subBand * TbaseA
    maxPossibleFreq = len(Freq);
    minFreq = int(i * subBand) #minFreq is integer minimum frequency, this subBand
    maxFreq = int((i + 1) * subBand) # maxFreq is integer maximum frequency, this subBand
    while  maxFreq < maxPossibleFreq:
        plot_num = str(i)

        #Now we create the filenames:

        filename = 'spec_%d.00_%d.00_%s_coherence_%s_and_%s' % (Freq[minFreq],Freq[maxFreq],startTime,CA,CB)

        #np.savetxt(path+filename+'.txt', Coh_Output[minFreq:maxFreq])
        np.savetxt(filename+'.txt', Coh_Output[minFreq:maxFreq], fmt=['%.6e','%.4f'])

        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.plot(Freq[minFreq:maxFreq:1], coh[minFreq:maxFreq:1], label = 'Coherence')
        ax1.set_xlabel('Frequency (Hz)',fontsize = 10)
        ax1.set_title(CA+' and '+CB+'; '+startTime+'; nAve= '+nAve, fontsize = 7)
        ax1.set_ylabel('Coherence',fontsize = 10)
        ax1.set_ylim([0,1])
        savefig(filename+'.png')
        #savefig(path+filename+'.png')
        clf()
        close

        # Setup to back to the top of the while loop.
        i += 1
        minFreq = int(i * subBand) #minFreq is integer minimum frequency, this subBand
        maxFreq = int((i + 1) * subBand) # maxFreq is integer maximum frequency, this subBand

    print('Done')
    # END def coherenceFromSFTs

# MAIN CODE STARTS HERE

#print sys.argv[0]

if len(sys.argv) < 3:
   print(' ')
   print('Find the coherence between SFTs in two specified directories ')
   print(' ')
   print('Usage: %s <pathToSFTsChanA> <pathToSFTsChanB> [subBand]' % sys.argv[0])
   print(' ')
   print('The optional subBand is the band in Hz to output in each plot. (Default is 100 Hz)')
   print(' ')
   exit(0)

#print sys.argv[1]
#print sys.argv[2]

pathToSFTsChannA = sys.argv[1]
pathToSFTsChannB = sys.argv[2]

if pathToSFTsChannA[-1] != '/':
   pathToSFTsChannA = pathToSFTsChannA + '/'

if pathToSFTsChannB[-1] != '/':
   pathToSFTsChannB = pathToSFTsChannB + '/'

# Default subBand
subBand=100

if len(sys.argv) >= 4:
   subBand =  sys.argv[3]

#coherenceFromSFTs('/home/kara.merfeld/public_html/fscan/test/output/fscans_2017_08_17_17_00_00_PDT_Thu/H1_GDS-CALIB_STRAIN/sfts/tmp/','/home/kara.merfeld/public_html/fscan/test/output/fscans_2017_08_17_17_00_00_PDT_Thu/H1_PEM-CS_ACC_BSC3_ITMX_X_DQ/sfts/tmp/')

coherenceFromSFTs(pathToSFTsChannA, pathToSFTsChannB, subBand)

