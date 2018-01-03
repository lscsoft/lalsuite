#!/usr/bin/python
import math, os, commands, shutil, sys, re
import matplotlib as matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_agg
import numpy as np
import argparse

# Plot time-frequency and second Fourier plane for TwoSpect
# Reads in tfdata.dat or ffdata.dat, made using TwoSpect --printData
# 02014-06-20 (JD 2456829)
# g m e a d o r s @ u m i c h .  e d u
# usage: ./createFFplane <pathtofile> <nSFTs> [optional: --fftf]

parser = argparse.ArgumentParser(description='Summarize the output files of a given Sco X-1 MDC data set')
parser.add_argument('path', type=str, help='Path to directory containing the output files')
parser.add_argument('Tobs', type=int, help='Observation time')
parser.add_argument('Tcoh', type=int, help='Coherence time, that is, one integration time of an SFT')
parser.add_argument('--fftf', action='store_true',help='Use to analyze frequency-frequency plane (default: time-frequency)')
args = parser.parse_args()

def tfplane(path, Tobs, Tcoh, fftf): 
        if fftf:
            print 'Printing data for frequency-frequency plane'
        else:
            print 'Printing data for time-frequency plane'

        # Read in data
        tfPlaneData = open(args.path,"r")
        tfList = []
        for ii, tfLine in enumerate(tfPlaneData):
            #if ii < 3*197:
            tfList.append(float(tfLine))
        tfArray = np.array(tfList)

        # Calculate expected number of SFTs
        nsft = int(np.floor(2 * Tobs / Tcoh)) - 1 
        print 'nsft = ' + str(nsft)
        # Calculate expected number of second Fourier transforms
        ntft = int(np.floor(nsft/2) + 1)
        print 'ntft = ' + str(ntft)

        # This really seems to be the right way to do it, with the transpose
        if fftf:
            tfShaped = np.reshape(tfArray, (len(tfArray)/ntft, ntft ))
        else:
            tfShaped = np.reshape(tfArray, (nsft, len(tfArray)/nsft )).T

        x, y = np.meshgrid(tfShaped[0, :], tfShaped[:, 0])
        if fftf:
            extensions = [0, ntft, 0, tfShaped.shape[0]]
        else:
            extensions = [0, nsft, 0, tfShaped.shape[0]]

        # Plot time-frequency plane
        fig = plt.figure()
        ax = fig.add_subplot(111)
        paramSpacePixelMap = ax.imshow(tfShaped, origin = 'lower', \
        interpolation = 'nearest', extent = extensions, cmap = 'jet')
        paramSpacePixelMap = fig.colorbar(paramSpacePixelMap, shrink = 0.5, extend = 'both')
        ax.set_aspect('auto')
        if fftf:
            ax.set_xlabel('2nd Frequency: f-prime (Hz)')
            ax.set_ylabel('Frequency bin: f (Hz) * Tcoh (s)')
            ax.set_title(\
            'Power in f-f plane ' + '\n' + 'Number of bins in data arrays (n f-prime, n f): ' +\
            str(tfShaped.T.shape) + ' \n '\
            ) 
        else:
            ax.set_xlabel('Time: SFT number (n)')
            ax.set_ylabel('Frequency bin: f (Hz) * Tcoh (s)')
            ax.set_title(\
            'Power in t-f plane ' + '\n' + 'Number of bins in data arrays (n t, n f): ' +\
            str(tfShaped.T.shape) + ' \n '\
            ) 
        if fftf:
            plt.savefig('ffplane.png')
            plt.savefig('ffplane.pdf')
        else:
            plt.savefig('tfplane.png')
            plt.savefig('tfplane.pdf')
        plt.close()
        plt.clf()

tfplane(args.path, args.Tobs, args.Tcoh, args.fftf)
