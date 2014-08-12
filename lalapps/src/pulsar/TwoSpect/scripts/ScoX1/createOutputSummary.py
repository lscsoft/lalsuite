#!/usr/bin/python
#from pywcsgrid2.allsky_axes import make_allsky_axes_from_header, allsky_header
#import pywcsgrid2.healpix_helper as healpix_helper
#from mpl_toolkits.basemap import Basemap
import math, os, commands, shutil, sys, re
import matplotlib as matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.backends.backend_agg
import numpy as np
import argparse

# Summarize the output files of a given Sco X-1 MDC data set
# 02014 -02-06 (JD 2456695)
# g m e a d o r s @ u m i c h .  e d u
# usage: ./createOutputSummary str(mdcVersion) str(observatory) str(pulsar).zfill(3)

parser = argparse.ArgumentParser(description='Summarize the output files of a given Sco X-1 MDC data set')
parser.add_argument('mdcVersion', type=str, help='Typically, "6"; alternately, for real data, use Science Run number')
parser.add_argument('observatory',type=str, help='Either H1, L1 or V1')
parser.add_argument('pulsar',type=str, help='Number of the pulsar in MDCv6; alternately, for real data, use floor(frequency in Hertz)')
parser.add_argument('--skyGrid',action='store_true', help='Produce sky maps for right ascension and declination')
parser.add_argument('--bypassSummary',action='store_true', help='Avoid regenerating the summary file, useful if already made')
parser.add_argument('--band',action='store_true', help='Use with 1 Hz templateSearch band, or standalone: make plots appropriate for a 5 Hz band')
parser.add_argument('--elsewhere', help='Path to directory containing the output files')
parser.add_argument('--massiveSummary',action='store_true', help='Deal with large output directories')
parser.add_argument('--plotSkyContour',type=float, help='Plot a circle around this point in the sky: two arguments, alpha and delta in radians', nargs=2)
parser.add_argument('--noiseTest',action='store_true', help='Make histrogram noise plots')
parser.add_argument('--templateSearch',action='store_true', help='Change commands in way needed for templateSearch')
parser.add_argument('--multiTemplateSearch',type=str, help='use instead of --templateSearch, specify number of bands in output directory')
parser.add_argument('--closed',action='store_true', help='Can be used, especially with --multiTemplateSearch, for working with closed pulsars')
parser.add_argument('--J1751',action='store_true',help='Use for a search on J1751 in real data')
parser.add_argument('--ScoX1S6',action='store_true',help='Use for a search on ScoX1 in S6 real data')
args = parser.parse_args()

def summarizer(mdcVersion, observatory, pulsar, args):
    headJobName = mdcVersion + '_' + observatory
    outdirectory = 'output_' + headJobName + "_pulsar-" + pulsar
    if args.elsewhere:
        outdirectory = args.elsewhere + outdirectory
    print outdirectory
    verboseSummaryFile = 'verbose_summary_' + headJobName + '-' + pulsar +'.txt'

    if args.bypassSummary:
        print 'Bypassing summary file creation'
    elif args.massiveSummary:
        print 'Taking alternate approach to a large output directory' 
        dagFileCore = headJobName
        if args.closed:
            dagFileCore = dagFileCore + '_closed'
        if args.elsewhere:
            jobDagFile = open(args.elsewhere + 'ScoX1_' + dagFileCore + '.dag', "r") 
        else:
            jobDagFile = open('ScoX1_' + dagFileCore + '.dag', "r") 
        outfilenameList = [] 
        for jobDagLine in jobDagFile:
            outfileLine = re.search("--outfilename=out_" + mdcVersion + \
            "_" + observatory + "_pulsar-" + pulsar + "_" + \
            "(?P<OUTFINT>\d+)\.(?P<OUTFFP>\d+)" + "_" + \
            "(?P<OUTDFINT>\d+)\.(?P<OUTDFFP>\d+)"+ ".dat", jobDagLine)
            if outfileLine:
                wholeOutfile = "out_" + mdcVersion + \
                "_" + observatory + "_pulsar-" + pulsar + "_" + \
                str(outfileLine.group(1) + '.' + outfileLine.group(2)) + '_' + \
                str(outfileLine.group(3) + '.' + outfileLine.group(4)) + '.dat'
                outfilenameList.append(wholeOutfile) 
        jobDagFile.close
        for ll, outfileEntry in enumerate(outfilenameList):
            if ll < 1e6: 
                if ll % 1e3 == 0:
                    print ll
                fileLocation = outdirectory + '/' + outfileEntry 
                grepCommand = 'grep -i ' + fileLocation + ' -e h0'
                os.system(grepCommand + ' >> ' + verboseSummaryFile)
        print 'Done'
    elif args.J1751:
        print 'Analyzing J1751'
        outJ1751directory = args.elsewhere + '/output_' + mdcVersion + \
        '_J1751_' + observatory
        outJ1751file = 'out_' + mdcVersion + '_J1751_' + observatory + \
        '_.dat'
        os.system('cat ' + outJ1751directory + '/' + outJ1751file + \
        ' | grep "h0 =" >> '+ verboseSummaryFile)
    elif args.ScoX1S6:
        print 'Analyzing ScoX1'
        outScoX1directory = args.elsewhere + '/output_S6' + \
        '_ScoX1_' + observatory + '-band-' + str(args.pulsar).zfill(4)
        outScoX1file = 'out_S6_ScoX1_' + observatory + \
        '-band-' + str(args.pulsar).zfill(4) + '*.dat'
        os.system('cat ' + outScoX1directory + '/' + outScoX1file + \
        ' | grep "h0 =" >> '+ verboseSummaryFile)
    else:
        # When already done, one need not repeat this step
        if args.skyGrid:
            print 'Producing output summary'
            listingOfSkyFiles = os.listdir(outdirectory)
            listingOfErrorFiles = [s for s in listingOfSkyFiles if 'err' in s]
            rangeOfErrorFiles = range(1, len(listingOfErrorFiles) + 1)
            os.system('rm ' + verboseSummaryFile)
            for iii in rangeOfErrorFiles:
                os.system('cat ' + outdirectory + '/TwoSpect.err.TwoSpect_' + str(iii) \
                + ' | grep "h0 =" >> '\
                + verboseSummaryFile)
        else:
            os.system('cat ' + outdirectory + '/out_' + headJobName + '* | grep "h0 =" > '\
            + verboseSummaryFile)
    
    # Now we have a single file containing all the output from all the
    # sky points around an injection. Let us read in the data
    verboseData = open(verboseSummaryFile, "r")
    rightAscensionList = []
    declinationList = []
    fList = []
    dfList = []
    RList = []
    ProbList = []
    for verboseLine in verboseData:
        verboseString = str(verboseLine)
        if args.skyGrid:
            # To save space, I use alpha to stand for right ascension
            # and delta for declination inside the regular expressions
            alphaLine = re.search("RA = (?P<ALPHASIGN>\-?)" \
            + "(?P<ALPHAINT>\d+)\.(?P<ALPHAFP>\d+)", \
            verboseString)
            rightAscensionList.append(float(alphaLine.group(1) + alphaLine.group(2) + \
            '.' + alphaLine.group(3)))
            deltaLine = re.search("DEC = (?P<DELTASIGN>\-?)" \
            + "(?P<DELTAINT>\d+)\.(?P<DELTAFP>\d+)", \
            verboseString)
            declinationList.append(float(deltaLine.group(1) + deltaLine.group(2) + \
            '.' + deltaLine.group(3))) 
        fLine = re.search("fsig = (?P<FSIGN>\-?)" \
        + "(?P<FINT>\d+)\.(?P<FFP>\d+)", \
        verboseString)
        fList.append(float(fLine.group(1) + fLine.group(2) + \
        '.' + fLine.group(3)))
        dfLine = re.search("df = (?P<DFSIGN>\-?)" \
        + "(?P<DFINT>\d+)\.(?P<DFFP>\d+)", \
        verboseString)
        dfList.append(float(dfLine.group(1) + dfLine.group(2) + \
        '.' + dfLine.group(3)))
        RLine = re.search("R = (?P<RSIGN>\-?)" \
        + "(?P<RINT>\d+)\.(?P<RFP>\d+)", \
        verboseString)
        RList.append(float(RLine.group(1) + RLine.group(2) + \
        '.' + RLine.group(3)))
        ProbLine = re.search("Prob = (?P<PROBSIGN>\-?)" \
        + "(?P<PROBINT>\d+)\.(?P<PROBFP>\d+)", \
        verboseString)
        ProbList.append(float(ProbLine.group(1) + ProbLine.group(2) + \
        '.' + ProbLine.group(3)))
    verboseData.close
 
    if args.skyGrid:
        # We should now have four lists, with as many
        # entries as there are sky points. Convert to NumPy arrays!
        rightAscensionArray = np.array(rightAscensionList)
        declinationArray = np.array(declinationList)
    fArray = np.array(fList)
   

    if (args.band or args.noiseTest) or \
    (args.templateSearch or args.multiTemplateSearch) or \
    args.J1751 or args.ScoX1S6:
        print 'Frequency list encompasses this many bins: ' + str(len(fArray))
    else:
        # Adjust f Array to represent frequency offset from injection
        fArray = fArray - fArray[(len(fArray)-1)/2]
    dfArray = np.array(dfList)
    if (args.band or args.noiseTest) or \
    (args.templateSearch or args.multiTemplateSearch) or \
    args.J1751 or args.ScoX1S6:
        print 'Frequency modulation list encompasses this many bins: ' + str(len(dfArray))
    RArray = np.array(RList)
    ProbArray = np.array(ProbList)

    if args.skyGrid:
        # By design, these are one-dimensional arrays that increment in
        # declination, then in right ascension, e.g., the (RA, dec) values go
        # (0, 0), (0, 1), (0, 2), (0, 3),... (1, 0), (1, 1), (1, 2), (1, 3)...
        # For convenience, we can reshape these arrays so we can use image plotters
        # We want a map where RA increases left to right and dec from bottom to top
        # i.e., where 
        # (0, 3), (1, 3), (2, 3), (3,3)
        # ...
        # (0, 0), (1, 0), (2, 0), (3,0)
        # To do this we first reshape by the length of right ascension, 
        # to break up the arrays in a matrix,
        # then transpose, to ensure RA increases the right way, 
        raLen = len(np.unique(rightAscensionArray))
        decLen = len(np.unique(declinationArray))
        raShaped = np.reshape(rightAscensionArray, (raLen, decLen)).T
        decShaped = np.reshape(declinationArray, (raLen, decLen)).T 
        # Having checked that this plots correctly (verified by substituting
        # raShaped or decShaped into the final plot), we then define the extents
        # of the plot so we have plot axis labels

        # We discover, however, that this method gets negative values of
        # declination backwards. That is tolerable for a local patch but
        # not for an all-sky map. How do we fix it? Sorting based on unique
        # values of declination is conceptually simplest:
        # This provides a key that we can use for all the matrices
        # Finally, it has to be flipped up-down
        print decShaped.shape
        I = np.argsort(decShaped[:,0])
        decShaped = decShaped[I,:]
        decShaped = decShaped[::-1,:]
        raShaped = raShaped[I,:]
        raShaped = raShaped[::-1,:]
    
        # All this shuffling is, of course, confirmed by plotting
        # ra and dec instead of Prob and R and checking that
        # the heat map plots make sense.
        x, y = np.meshgrid(raShaped[0, :], decShaped[:, 0])
        extensions = [x[0, 0], x[-1, -1], y[-1, -1], y[0, 0]]
        #extensions = [0, 2*np.pi, -np.pi/2, np.pi/2]
        #print extensions


        # Now we simply plot the probability and R statistic

        ProbShaped = np.absolute(np.reshape(ProbArray, (raLen, decLen)).T)
        ProbShaped = ProbShaped[I,:]
        ProbShaped = ProbShaped[::-1,:]
        RShaped = np.reshape(RArray, (raLen, decLen)).T
        RShaped = RShaped[I,:]
        RShaped = RShaped[::-1,:]

        if args.plotSkyContour:
            # Thanks to http://web.archiveorange.com/archive/v/74jMQNS3vyH1xzrlqv6c
            # for tips on how to plot an all-sky map
            # The formula for the great-circle distance between points on a sphere is,
            # per the Wikipedia article on Great-circle distance,
            # d = r * dsigma; we take r = 1 and dsigma as below
            # dsigma = arccos(sin(phi_s)*sin(phi_f) + cos(phi_s)*cos(phi_f)*cos(abs(lambda_s-lambda_f)))
            # where the phi denotes a latitude, lambda a longitude, and s and f are the two points
            def sphereDistance(phi_s, lambda_s, phi_f, lambda_f):
                return np.arccos(np.sin(phi_s) * np.sin(phi_f) + \
                np.cos(phi_s) * np.cos(phi_f) * np.cos(np.absolute(lambda_s - lambda_f)))
            # For MDCv2 vers MDCv5 study:
            # Earth direction approximately toward RA 135 degrees, lat plus 16, based
            # on solar direction three-months prior (did not use full JPL correction)
            gcDistanceTrue = []
            gcDistance = [sphereDistance(args.plotSkyContour[1], args.plotSkyContour[0], lat, lon) for lat in y[:,0] for lon in x[0,:]]
        # Plot the sky map for log10probability
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        m = Basemap(projection='moll',lon_0=180, lat_0=0)
        mapx, mapy = m(180/np.pi*x, 180/np.pi*y)
        skyMap = m.pcolormesh(mapx, mapy, ProbShaped, \
        ) #clip_on=False, interpolation='nearest', cmap='jet', extent=extensions,aspect=0.5)
        ax.grid('True')
        m.drawmeridians(np.arange(0,360,30), labels=[1,0,0,0],fontsize=10)
        m.drawparallels(np.arange(-90,90,30),  fontsize=10)
        skyMap = fig.colorbar(skyMap, extend = 'both', shrink=0.3)
        ax.set_xlabel('Right ascension (rad)')
        ax.set_ylabel('Declination (rad)')
        if args.plotSkyContour:
            contourBlurb = '\n Great circle distance contours centered around'\
            + '\n (alpha, delta)= ' + \
            '(' + str(args.plotSkyContour[0]) + ', ' + \
            str(args.plotSkyContour[1]) + ')'
        else:
            contourBlurb = ''
        ax.set_title('All-sky map of TwoSpect' + \
        '\n log 10 probability vs sky position for ' + headJobName + \
        '\n Note: RA is 0h on left, 12h center, increasing left-to-right' + \
        contourBlurb)
        # The following title was only used for the MDCv2 vs MDCv5 study:
        #ax.set_title('all-sky map of TwoSpect' + \
        #'\n log 10 probability vs sky position for ' + headJobName + \
        #'\n with constant-Doppler-correction contours' + \
        #'\n approximating v_Earth as toward solar direction' +\
        #'\n three months prior to GPS time 815432000 (2005-11-07), i.e.,'
        #'\n RA 9h, dec +16 deg'
        #'\n (note: RA is 0h on left, 12h center, increasing left-to-right)')
        if args.plotSkyContour:
            print 'Plotting sky contour around this (alpha, delta) location: ' + \
            str(args.plotSkyContour[0]) + ' ' + \
            str(args.plotSkyContour[1])
            gcDistanceArray = np.asarray(gcDistance).reshape(y.shape[1],-1)
            m.contour(mapx,mapy,gcDistanceArray)
        canvas = matplotlib.backends.backend_agg.FigureCanvasAgg(fig)
        canvas.print_figure("all-sky-map_prob_" + headJobName + ".png")
        canvas.print_figure("all-sky-map_prob_" + headJobName + ".pdf")
        plt.close()

        #Plot the sky map for R
        fig = plt.figure(figsize=(8,8))
        ax = fig.add_axes([0.1,0.1,0.8,0.8])
        m = Basemap(projection='moll',lon_0=180, lat_0=0)
        #ax = fig.add_subplot(1, 1, 1, projection='mollweide')
        mapx, mapy = m(180/np.pi*x, 180/np.pi*y)
        skyMapR = m.pcolormesh(mapx, mapy, RShaped, \
        ) # interpolation = 'nearest', cmap='jet', extent=extensions)
        skyMapR = fig.colorbar(skyMapR, extend = 'both')
        ax.grid('True')
        #ax.set_aspect('auto')
        #ax.set_xlim((0,2*np.pi))
        ax.set_xlabel('Right ascension (rad)')
        ax.set_ylabel('Declination (rad)')
        ax.set_title('R statistic vs sky position for ' + headJobName)
        fig.savefig('all-sky-map_R_' + headJobName + '.png')
        fig.savefig('all-sky-map_R_' + headJobName + '.pdf')
        plt.close()
    else:
        # The above method should also work for f and df,
        # with f increasing horizontally, df vertically,
        # as in figure 4 of the TwoSpect methods paper
        fLen = len(np.unique(fArray))
        dfLen = len(np.unique(dfArray))
        if (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            # This override is necessary when using templateSearch because the
            # df are, for it, all unique. However, this makes it tricky because
            # the graph really is skewed
            dfLen = len(fArray)/fLen
            if args.multiTemplateSearch:
                fLenOne = fLen/int(args.multiTemplateSearch)
                dfLenCml = range(0, int(args.multiTemplateSearch))
                dfLenPer = range(0, int(args.multiTemplateSearch))
                startCount = range(0, int(args.multiTemplateSearch))
                midCount = range(0, int(args.multiTemplateSearch))
                fStart = float(fArray[0])
                fEnd = float(fArray[-1])
                fSpan = (float(fEnd)-float(fStart))/float(args.multiTemplateSearch)
                print 'fSpan (Hz) of each band is ' + str(fSpan)
                print 'Starting frequency (Hz) is ' + str(fStart)
                print 'Ending frequency (Hz) is ' + str(fEnd)
                for nn in range(0, int(args.multiTemplateSearch)):
                   counterLeft = int( (int(nn) + 1)/float(fSpan))
                   counterRight = int( int(args.multiTemplateSearch)/float(fSpan))
                   fpTolerance = 2e-13
                   if counterLeft < counterRight:
                       dfLenCml[nn] = len([ff for ff in fArray if (ff < float(fStart)+float(nn+1)*fSpan)])
                       #dfLenPer[nn] = len([ff for ff in fArray if (ff >= fStart + nn) and (ff < fStart+nn+1)])
                       startCount[nn] = len([ff for ff in fArray if abs(ff - (float(fStart) + float(nn)*fSpan) ) < fpTolerance])
                       midCount[nn] = len([ff for ff in fArray if abs(ff - (float(fStart) + float(nn+0.5)*fSpan) ) < fpTolerance])
                   elif counterLeft == counterRight:
                       print 'Ending after this many Hertz ' + str((nn+1)*fSpan)
                       dfLenCml[nn] = len([ff for ff in fArray if (ff <= float(fStart)+float(nn+1)*fSpan)])
                       #dfLenPer[nn] = len([ff for ff in fArray if (ff >= fStart + nn) and (ff <= fStart+nn+1)])
                       startCount[nn] = len([ff for ff in fArray if abs(ff - (float(fStart) + float(nn)*fSpan) ) < fpTolerance])
                       midCount[nn] = len([ff for ff in fArray if abs(ff - (float(fStart) + float(nn+0.5)*fSpan) ) < fpTolerance])
                   else:
                       print 'Malfunction, reached nn of ' + str(nn)
                print 'startCount, midCount, dfLenCml:'
                print startCount
                print midCount
                print dfLenCml
                for qq in range(0, int(args.multiTemplateSearch)):
                    if qq == 0:
                        startIndex = 0
                        endIndex = int(dfLenCml[qq])
                    else:
                        startIndex = int(dfLenCml[qq-1]) + startCount[qq] - midCount[qq]
                        endIndex = int(dfLenCml[qq])
                    if qq + 1 == int(args.multiTemplateSearch):
                        fLenSub = fLenOne + 1
                    else:
                        fLenSub = fLenOne
                    # Special contigency; sometimes we want the same code
                    # even though multiTemplateSearch is not really needed
                    if int(args.multiTemplateSearch) == 1:
                        fLenSub = fLenOne
                    #print fLenSub
                    #print midCount[qq]
                    #print str(fArray[startIndex])
                    #print str(fArray[endIndex-1])
                    fSubShaped = np.reshape(fArray[startIndex:endIndex], (fLenSub, midCount[qq])).T
                    dfSubShaped = np.reshape(dfArray[startIndex:endIndex], (fLenSub, midCount[qq])).T
                    ProbSubShaped = np.absolute(np.reshape(ProbArray[startIndex:endIndex], (fLenSub, midCount[qq])).T)
                    RSubShaped = np.reshape(RArray[startIndex:endIndex], (fLenSub, midCount[qq])).T
                    #print dfSubShaped.shape
                    noRowToAdd = np.max(midCount) - dfSubShaped.shape[0]
                    fRowToAdd = fSubShaped[0,:]
                    dfRowToAdd = dfSubShaped[-1,:]
                    if noRowToAdd > 0:
                        #print (np.zeros((noRowToAdd, fLenSub))).shape
                        for qqq in range(0, noRowToAdd):
                            fSubShaped = np.vstack((fSubShaped, fRowToAdd))
                            dfSubShaped = np.vstack((dfSubShaped, dfRowToAdd))
                            #ProbSubShaped = np.vstack((ProbSubShaped, fRowToAdd))
                            #RSubShaped = np.vstack((RSubShaped, dfRowToAdd))
                        ProbSubShaped = np.vstack((ProbSubShaped, np.zeros((noRowToAdd, fLenSub))))
                        RSubShaped = np.vstack((RSubShaped, np.zeros((noRowToAdd, fLenSub))))
                    if qq == 0:
                        fShaped = fSubShaped
                        dfShaped = dfSubShaped
                        ProbShaped = ProbSubShaped
                        RShaped = RSubShaped
                    else:
                        fShaped = np.hstack((fShaped, fSubShaped))
                        dfShaped = np.hstack((dfShaped, dfSubShaped))
                        ProbShaped = np.hstack((ProbShaped, ProbSubShaped))
                        RShaped = np.hstack((RShaped, RSubShaped))
                    #print dfShaped.shape
            if args.J1751 or \
            args.noiseTest:
                fShaped = np.reshape(fArray, (fLen, dfLen)).T
                dfShaped = np.reshape(dfArray, (fLen, dfLen)).T
            elif args.ScoX1S6:
                print 'Shape of real-data Scorpius X-1 array: '
                print str(dfShaped.shape)
        else:
            fShaped = np.reshape(fArray, (fLen, dfLen)).T
            dfShaped = np.reshape(dfArray, (fLen, dfLen)).T
        if (args.band or args.noiseTest) or \
        (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6: 
            print 'Number of bins in data arrays: ' + str(fShaped.shape)
        x, y = np.meshgrid(fShaped[0, :], dfShaped[:, 0])
        extensions = [x[0, 0], x[-1, -1], y[0, 0], y[-1, -1]]
        if not args.multiTemplateSearch:
            ProbShaped = np.absolute(np.reshape(ProbArray, (fLen, dfLen)).T)
            RShaped = np.reshape(RArray, (fLen, dfLen)).T
      
        if (args.band or args.noiseTest) or \
        (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            ProbCenter = ProbShaped.max()
            RCenter = RShaped.max()
            centerString = 'maximum value: '
        else:
            ProbCenter = ProbShaped[(np.shape(ProbShaped)[0]-1)/2, (np.shape(ProbShaped)[1]-1)/2]
            RCenter = RShaped[(np.shape(RShaped)[0]-1)/2, (np.shape(RShaped)[1]-1)/2]
            centerString = 'center value: '
        centerProbSpotF = str(fShaped.compress((ProbShaped == ProbCenter).flat)[0])
        centerProbSpotDF = str(dfShaped.compress((ProbShaped == ProbCenter).flat)[0])
        centerRSpotF = str(fShaped.compress((RShaped == RCenter).flat)[0])
        centerRSpotDF = str(dfShaped.compress((RShaped == RCenter).flat)[0])
        # Plot probability two-dimensionally for f and df
        if (args.band or args.noiseTest) or \
        (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            fig = plt.figure(figsize=(12,12))
        else:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        paramSpacePixelMap = ax.imshow(ProbShaped, origin='lower', \
        interpolation = 'nearest', extent = extensions, cmap = 'jet')
        paramSpacePixelMap = fig.colorbar(paramSpacePixelMap, shrink = 0.5, extend = 'both')
        if (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            print 'Skipping probability grid'
        else:
            ax.grid('True')
        if (args.band or args.noiseTest) or \
        (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            ax.set_aspect('auto')
            ax.set_xlabel('Frequency: f (Hz)')
        else:
            ax.set_xlabel('Frequency offset from injection: f (Hz)')
        ax.set_ylabel('Modulation depth: df (Hz)')
        if args.ScoX1S6:
            ax.set_title('False alarm probability, absolute value of base 10 logarithm: \n \
            log10p vs parameters for band starting ' + pulsar + ' Hz at ' + observatory + ' \n \
            ' + centerString + str(ProbCenter) + ' at (df, f) = (' + centerProbSpotDF +', ' + centerProbSpotF + ') Hz \n \
            Number of bins in data arrays (df, f): ' + str(fShaped.shape) + ' \n \
            ')
        else:
            ax.set_title('False alarm probability, absolute value of base 10 logarithm: \n \
            log10p vs parameters for pulsar ' + pulsar + ' at ' + observatory + ' \n \
            ' + centerString + str(ProbCenter) + ' at (df, f) = (' + centerProbSpotDF +', ' + centerProbSpotF + ') Hz \n \
            Number of bins in data arrays (df, f): ' + str(fShaped.shape) + ' \n \
            ')
        plt.savefig('DFvsFresultsProb-' + observatory + '_pulsar-' + pulsar + '.png')
        plt.savefig('DFvsFresultsProb-' + observatory + '_pulsar-' + pulsar + '.pdf')
        plt.close()
        plt.clf()

        # Plot R statistic two-dimensionally for f and df
        if (args.band or args.noiseTest) or \
        (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            fig = plt.figure(figsize=(12,12))
        else:
            fig = plt.figure()
        ax = fig.add_subplot(111)
        paramSpacePixelMap = ax.imshow(RShaped, origin = 'lower', \
        interpolation = 'nearest', extent = extensions, cmap = 'jet')
        paramSpacePixelMap = fig.colorbar(paramSpacePixelMap, shrink = 0.5, extend = 'both')
        if (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            print 'Skipping R grid'
        else:
            ax.grid('True')
        if (args.band or args.noiseTest) or \
        (args.templateSearch or args.multiTemplateSearch) or \
        args.J1751 or args.ScoX1S6:
            ax.set_aspect('auto')
            ax.set_xlabel('Frequency: f (Hz)')
        else:
            ax.set_xlabel('Frequency offset from injection: f (Hz)')
        ax.set_ylabel('Modulation depth: df (Hz)')
        ax.set_title('R statistic \n \
        R vs parameters for pulsar ' + pulsar + ' at ' + observatory + ' \n \
        ' + centerString + str(RCenter) + ' at (df, f) = (' + centerRSpotDF +', ' + centerRSpotF + ') Hz \n \
        Number of bins in data arrays (df, f): ' + str(fShaped.shape) + ' \n \
        ') 
        plt.savefig('DFvsFresultsR-' + observatory + '_pulsar-' + pulsar + '.png')
        plt.savefig('DFvsFresultsR-' + observatory + '_pulsar-' + pulsar + '.pdf')
        plt.close()
        plt.clf()

if args.J1751 or args.ScoX1S6:
  summarizer('S' + args.mdcVersion, args.observatory, str(args.pulsar).zfill(3), args)
else:
  summarizer('mdcv' + args.mdcVersion, args.observatory, str(args.pulsar).zfill(3), args)
