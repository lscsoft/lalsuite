#!/usr/bin/python
import os, commands, shutil, sys, re
import numpy as np
import argparse

# Create a Condor DAG submit file to analyze sky locations and parameter space
# in the Sco X-1 TwoSpect mock data challenge
# 02013-11-18 (JD 2456614)
# g m e a d o r s @ u m i c h .  e d u
# usage: ./createScoX1dag.py str(observatory) int(fSteps) int(dfSteps)

parser = argparse.ArgumentParser(description='Create a Condor DAG submit file to analyze sky locations')
parser.add_argument('observatory', type=str, help='Either H1, L1 or V1')
parser.add_argument('fSteps', type=int, help='Number of frequency steps')
parser.add_argument('dfSteps', type=int, help='Number of frequency modulation steps')
parser.add_argument('--singleBand', type=int, help='Select a single 5 Hz band from the MDCv6 data set')
parser.add_argument('--skyGrid', type=int, help='Run over a grid of all sky locations')
parser.add_argument('--onlyOne', type=int, help='Choose only one pulsar from the list')
parser.add_argument('--noiseTest', type=float, help='Sample all bands for noise at this frequency resolution (Hz)')
parser.add_argument('--templateSearch', type=int, help='Use the templateSearch option for a 5 Hz band. Specify a pulsar number')
parser.add_argument('--templateSearchOpen', action='store_true', help='templateSearch, run all open pulsar bands, but not closed')
parser.add_argument('--templateSearchClosed', action='store_true', help='templateSearch, run all closed bands, but not open')
args = parser.parse_args()
if args.singleBand or args.templateSearch: 
    print 'Looking only in a 5 Hz band from the MDCv6 data set at the following pulsar: ' + str(args.singleBand)
elif args.templateSearchOpen or args.templateSearchClosed:
    print 'Looking in 5 Hz bands from the MDCv6 data set.'
if args.skyGrid:
    print 'Running over a grid of multiple sky locations'

def sftFileListParts():
    #For the closed search
    # They should correspond to the (floor of the five-Hertz band, minus
    # half a Hertz) times Tcoh
    sftFileListPart1Closed = [41580, 49980, 58380, 75180, 125580, 150780, 159180,\
    175980, 192780, 201180, 209580, 226380, 251580, 276780]
    sftFileListPart2Closed = [129420, 140220, 143820, 161820, 172620, 183420,\
    187020, 194220, 197820, 201420, 212220, 223020, 230220, 233820, 237420,\
    241020, 244620, 248220, 251820, 255420, 262620, 266220, 269820, 273420,\
    287820, 291420, 295020, 309420, 316620, 334620, 388620, 395820, 399420,\
    428220, 475020, 493020]

    # Manually construct the sftFileList, because it is so idiosyncratic 
    sftFileListPart1 = [66780, 83580, 91980, 100380, 108780, 117180, 133980,\
    142380, 167580, 184380, 217980, 234780, 243180, 259980, 268380, 285180,\
    293580]
    sftFileListPart2 = [133020, 136620, 147420, 151020, 154620, 158220,\
    165420, 169020, 176220, 179820, 190620, 208620, 215820, 219420, 277020,\
    298620, 320220, 331020, 345420, 374220, 377820, 403020, 417420, 421020,\
    424620, 449820, 460620, 464220, 467820, 478620, 489420, 500220, 521820]
    return [sftFileListPart1, sftFileListPart2, \
    sftFileListPart1Closed, sftFileListPart2Closed]


def sftNameMaker(observatory, tcoh, binname, args):
    if args.templateSearchClosed:
        openFlag = "closed"
    else:
        openFlag = "open"
    sftFileRoot = "/home/egoetz/TwoSpect/scox1_mdc6/sfts/" + \
    observatory + "/" 
    sftFilePart1 = "s_sfts/" + openFlag + "/" + observatory[0] + "-" + observatory[1] + \
    "_" + observatory + "_"
    sftFilePart2 = "SFT_SCO_X1_MDCv6_"
    # EXAMPLE: print "/home/egoetz/TwoSpect/scox1_mdc6/sfts/H1/840s_sfts/H-1_H1_840SFT_SCO_X1_MDCv6-313.5Hz_263340"
    return sftFileRoot + str(tcoh) + sftFilePart1 + str(tcoh) + sftFilePart2 + str(binname)
 
def tablereader(tableName, observatory, args):
    print 'Name of the MDC table:'
    print tableName
    tableData = open(tableName, "r")
        
    for k, tableLine in enumerate(tableData):
        if k == 0:
            if args.templateSearchClosed:
                # Note this would need to be moved outside if the 
                # first pulsar were not closed. In that case, just do
                # a check to see if these lists exist
                raInjList = np.asarray(4.275699238500, dtype=np.float64)
                decInjList = np.asarray(-0.272973858335, dtype=np.float64)
                fInjList = np.asarray(tableLine.split()[1], dtype=np.float64)
                pulsarNoList = np.asarray(tableLine.split()[0], dtype=np.int)
                if np.asarray(tableLine.split()[1], dtype=np.float32) < np.asarray(360.0):
                    TcohLine = 840
                else:
                    TcohLine = 360
                TcohList = (np.asarray(TcohLine)).astype(int)
                PList = [str(68023.7)]
                asiniList = [str(1.44)]
            else:
                raInjList = np.zeros((1,1))
                decInjList = np.zeros((1,1))
                fInjList = np.zeros((1,1))
                pulsarNoList = np.zeros((1,1))
                TcohList = np.zeros((1,1))
                PList = ['0']
                asiniList = ['0']
        if k > 0:
            if args.templateSearchClosed:
                if tableLine.split()[4] == 'closed':
                    raInjList = np.vstack([raInjList, np.asarray(4.275699238500, dtype=np.float64)])
                    decInjList = np.vstack([decInjList, np.asarray(-0.272973858335, dtype=np.float64)])
                    fInjList = np.vstack([fInjList, np.asarray(tableLine.split()[1], dtype=np.float64)])
                    pulsarNoList = np.vstack([pulsarNoList, np.asarray(tableLine.split()[0], dtype=np.int)])
                    if np.asarray(tableLine.split()[1], dtype=np.float32) < np.asarray(360.0):
                        TcohLine = 840
                    else:
                        TcohLine = 360
                    TcohList = np.vstack([TcohList, np.asarray(TcohLine)]).astype(int)
                    PList.append(str(68023.7))
                    asiniList.append(str(1.44))
            else:
                raInjList = np.vstack([raInjList, np.asarray(tableLine.split()[1], dtype=np.float64)])
                decInjList = np.vstack([decInjList, np.asarray(tableLine.split()[3], dtype=np.float64)])
                fInjList = np.vstack([fInjList, np.asarray(tableLine.split()[4], dtype=np.float64)])
                pulsarNoLine = re.search("Pulsar (?P<PULSARINT>\d+)", tableLine)
                pulsarNoList = np.vstack([pulsarNoList, np.asarray(np.asarray(pulsarNoLine.group(1)))])
                if np.asarray(tableLine.split()[4], dtype=np.float32) < np.asarray(360.0):
                    TcohLine = 840
                else:
                    TcohLine = 360
                TcohList = np.vstack([TcohList, np.asarray(TcohLine)]).astype(int) 
                PList.append(str(tableLine.split()[11]))
                asiniList.append(str(tableLine.split()[9]))

    sftFileList = []
    sftFileListPartList = sftFileListParts()
    if args.templateSearchClosed:
        sftFileListPart1 = sftFileListPartList[2]
        sftFileListPart2 = sftFileListPartList[3]
    else:
        sftFileListPart1 = sftFileListPartList[0]
        sftFileListPart2 = sftFileListPartList[1]
    [sftFileList.append(sftNameMaker(observatory, 840, binnumber, args))\
    for binnumber in sftFileListPart1]
    [sftFileList.append(sftNameMaker(observatory, 360, binnumber, args))\
    for binnumber in sftFileListPart2]
            
    if args.templateSearchClosed:
        print 'Done sorting closed pulsar table.'
    else:
        raInjList = raInjList[1:]
        decInjList = decInjList[1:]
        fInjList = fInjList[1:]
        pulsarNoList = pulsarNoList[1:]
        TcohList = TcohList[1:]
        PList = PList[1:]
        asiniList = asiniList[1:]
          
    return [raInjList, decInjList, fInjList, pulsarNoList, TcohList, sftFileList, PList, asiniList]
    


def categorizer(Tcoh, raInj, decInj, fInj, observatory, pulsarNo, sftFile, jobInc, Period, asini, fSteps, dfSteps, args):
    # Define a function to edit file objects conveniently
    def h(text):
        result = condorObject.write(text + '\n')
        return result
    # Define a second, similar function to write to a different file
    # Keep them separate for safety
    def g(text):
        result = dagObject.write(text + '\n')
        return result

    # Specify user name and the parent directory for analysis
    username = "grant.meadors"
    userDirectory = "/home/" + username + "/ScoX1/"
    parentDirectory = os.getcwd()

    # Specify which Sco X-1 MDC data set and observatory are to be used
    mdcVersion = "mdcv" + "6"
    headJobName = mdcVersion + '_' + observatory + '_' + pulsarNo
    headJobShort = mdcVersion + '_' + observatory
    if args.templateSearchClosed:
        headJobShort = headJobShort + "_closed"

    # Specify the right ascension and declination
    # where the signal is intended to be injected
    # For testing
    #raInj = 0.0
    #decInj = 0.0
    # MDC v2
    #raInj = 3.38
    #decInj = 0.19
    # MDC v3
    #raInj = 2.20 
    #decInj = -1.01
    # MDC v4
    #raInj = 3.42
    #decInj = 0.45
    # For all-sky testing
    #raInj = 3.14
    #decInj = 0.00
    # MDC v6
    #raInj = 4.275699238500	
    #decInj = -0.27
    # MDC v6 should take alpha and dec from input arguments
    # Using the argument parser for the all-sky test
    if args.skyGrid:
        raInj = 3.14
        decInj = 0.00

        # Specify how broadly to search around the injection point
        # The number of sky points will be (2*steps+1)^2
        interval = 1/8.0 * 25.0/args.skyGrid
        steps = args.skyGrid
        bound = steps * interval
        raRangeRaw = [-1.0* bound + interval * x for x in range(0, 2*steps +1)]
        #For the all-sky map, Dec will be half the spacing of RA
        decRangeRaw = [-0.5* bound + interval/2.0 * x for x in range(0, 2*steps +1)]
        raRange = [raInj + y for y in raRangeRaw]
        decRange = [decInj + y for y in decRangeRaw]
    else:
        raRange = raInj
        decRange = decInj
    if args.skyGrid:
        print 'Right ascension range (radians):'
        print raRange
        print 'Declination range (radians):'
        print decRange
    if (args.singleBand or args.templateSearch or \
    args.templateSearchOpen or args.templateSearchClosed) or args.noiseTest:
        sftFileListPartList = sftFileListParts()
        
        def FloorBinFinder(Tcoh, fInj, sftFile, sftFileSubList):
            fInjBin = Tcoh*fInj
            fBinMatchList = []
            [fBinMatchList.append(bin) for bin in sftFileSubList if bin <= fInjBin]
            # Error checking: the sftFile and the binList information match?
            if int(fBinMatchList[-1]) == int(sftFile.split('_')[-1]):
                 return fBinMatchList[-1]
            else:
                print 'Error: SFT bins for this pulsar appear incorrect'
                return None
        if Tcoh == 840:
            floorFlag = 0
        if Tcoh == 360:
            floorFlag = 1
        if args.templateSearchClosed:
            floorFlag = floorFlag + 2
        fFloor = FloorBinFinder(Tcoh, fInj, sftFile, \
        sftFileListPartList[floorFlag])    
        if args.noiseTest:
            fInterval =  args.noiseTest
        elif args.singleBand or args.templateSearch or \
        args.templateSearchOpen or args.templateSearchClosed:
            fInterval = 1 / (2 * float(Tcoh))
        else:
            print 'Conditional statements may be mixed-up'
        dfInterval = 1 / (4 * float (Tcoh))
        fBound = 5
        fSteps = int(fBound / fInterval)
        # Generate an odd number of steps
        if args.noiseTest:
            fRange = [0.5 + (fFloor/float(Tcoh)) + fInterval*y for y in range(0, fSteps + 1) ]
        elif args.templateSearch or \
        args.templateSearchOpen or args.templateSearchClosed:
            fRange = [fFloor/float(Tcoh)+ y for y in range(1,6)]
        else:
            fRange = [0.5 + (fFloor/float(Tcoh)) + fInterval*y for y in range(0, fSteps - 1) ]
        dfFloor = 2 * np.pi * (fFloor/Tcoh) * float(asini) / float(Period)
        dfCeiling = 2 * np.pi * (fRange[-1]) * float(asini) / float(Period)
        dfBound = dfCeiling - dfFloor
        dfAvg = (dfCeiling + dfFloor)/2
        dfSteps = int(dfBound / dfInterval)
        # Ensure that the modulation range search includes at least five points,
        # centered on the most likely modulation depth for the band
        if args.noiseTest:
            dfScope = range(0, 1)
        elif args.templateSearch or \
        args.templateSearchOpen or args.templateSearchClosed:
            dfScope = range(0, 1)
        else:
            dfScope = range(0, dfSteps + 5)
        dfRange = [dfAvg + dfInterval*(-0.5*dfScope[-1] + y) for y in dfScope]
    else:
        # Specify the range of frequency to search
        fHypothesis = fInj
        fInterval = 1 / (2 * float(Tcoh))
        fBound = (fSteps - 1) * fInterval
        fRange = [fHypothesis - 0.5*fBound + fInterval*y for y in range(0, fSteps)]
        # Specify the range of frequency modulation to search
        dfHypothesis = 2 * np.pi * fInj  * float(asini) / float(Period)
        dfInterval = 1 / (4 * float (Tcoh))
        dfBound = (dfSteps - 1) * dfInterval 
        dfRange = [dfHypothesis -0.5*dfBound + dfInterval*y for y in range(0, dfSteps)]

    # Choose a TwoSpect version
    # Grant David Meadors's latest development version
    TwoSpectVersion = \
    "/home/" + username + "/TwoSpect/dev/bin/lalapps_TwoSpect"
    # Grant David Meadors's version 1.1.27
    #TwoSpectVersion = \
    #"/home/" + username + "/master/opt/lscsoft/lalapps/bin/lalapps_TwoSpect"
    # Evan Goetz's version 1.1.27
    #TwoSpectVersion = \
    #"/home/egoetz/opt/lscsoft/bin/lalapps_TwoSpect"

    # Make a directory for the output logs
    os.system('mkdir -p output_' + headJobName)

    # Write a Condor sub file
    condorObject = open(parentDirectory + "/ScoX1_" + \
    headJobName + ".sub", "a")

    # Insert the contents of the file
    h("universe = vanilla")
    h("executable = " + TwoSpectVersion)
    h("output = " + "output_" + headJobName + "/TwoSpect.out.$(tagstring)")
    h("error = output_" + headJobName + "/TwoSpect.err.$(tagstring)")
    h("log = output_" + headJobName + "/TwoSpect.dag.log")
    h("request_memory = 3 GB")
    h("notification = never")
    h("environment = HOME=/home/" + username)
    h("")
    h("arguments = $(argList)")
    h("queue 1")
    h("")

    # Close the Condor sub file
    condorObject.close

    # Write each job of the DAG file
    # in the dagWriter function
    dagObject = open(parentDirectory + "/ScoX1_" + \
    headJobShort + ".dag", "a")

    if args.skyGrid:
        for ra in raRange:
            for dec in decRange:
                [dagWriter(g, observatory, headJobName, jobInc + 1 + dfRange.index(df) + len(dfRange)*(fRange.index(f) + len(fRange)*(decRange.index(dec) +len(decRange)*raRange.index(ra))), ra, dec, f, df, Tcoh, sftFile, headJobShort, Period) for f in fRange for df in dfRange]
    else:
        [dagWriter(g, observatory, headJobName, jobInc + 1 + dfRange.index(df) + len(dfRange)*fRange.index(f), raRange, decRange, f, df, Tcoh, sftFile, headJobShort, Period) for f in fRange for df in dfRange]
    dagObject.close

def dagWriter(g, observatory, headJobName, jobNumber, rightAscension, declination, f, df, Tcoh, sftFile, headJobShort, Period):
    # subtract half the fspan, 0.125 Hz, and round fmin
    # to the nearest 0.125 Hz
    # to supply the fspan argument

    #startTime is conditional on the observatory
   
    if str(Tcoh) == str(360):
        if str(observatory) == 'H1':
            startTime = str(1230338520)
        elif str(observatory) == 'L1':
            startTime = str(1230342660)
        elif str(observatory) == 'V1':
            startTime = str(1230337080)
        else:
            print 'Observatory choice (argument 1) not understood: choose H1, L1, or V1'
    elif str(Tcoh) == '840':
        if str(observatory) == 'H1':
            startTime = str(1230338760)
        elif str(observatory) == 'L1':
            startTime = str(1230342540)
        elif str(observatory) == 'V1':
            startTime = str(1230337080)
        else:
            print 'Observatory choice (argument 1) not understood: choose H1, L1, or V1'
    else:
        print 'Coherence time incorrectly assigned: must be 360 or 840 s'
        
    if args.templateSearch or \
    args.templateSearchOpen or args.templateSearchClosed:
        fminString = str(1*np.math.floor(f/1))
    else:
        fminString = str(np.math.floor(8*(f-0.125))/8)     
    if args.templateSearch or \
    args.templateSearchOpen or args.templateSearchClosed:
        #templateStringSet = ' --templateSearch' + \
        #' --fspan=' + str(1 - 1/(2*float(Tcoh)))
        templateStringSet = ' --templateSearch' + \
        ' --fspan=' + str(1)
    else:
        templateStringSet = ' --templateTest' + \
        " --templateTestF=" + str(f) + \
        " --templateTestDf=" + str(df) + \
        " --templateTestP=" + str(Period) + \
        ' --fspan=0.125'

    argumentList = \
    '"' + \
    ' --config=config_file_mdcv6.txt' + \
    " --t0=" + str(startTime) + \
    " --fmin=" + fminString + \
    " --skyRegion=(" + str(rightAscension)+ ',' + str(declination) + ")" + \
    templateStringSet + \
    " --Tcoh=" + str(Tcoh) + \
    " --SFToverlap=" + str(Tcoh/2) + \
    " --sftFile=" + str(sftFile) + \
    " --IFO=" + str(observatory) + \
    " --outfilename=out_" + headJobName + '_' + \
    str(f) + '_' + \
    str(df) + \
    ".dat"  + " --outdirectory=output_" + headJobName + '"' 
    tagStringLine = "TwoSpect_" + str(jobNumber)
    g("JOB " + tagStringLine + " ScoX1_" + headJobName + ".sub")
    g("VARS " + tagStringLine + " argList=" + argumentList + " tagString=" + '"' + tagStringLine + '"')

observatoryChoice = args.observatory

if args.templateSearchClosed:
    print 'Preparing MDCv6 closed search, reading band table...'
    tableName = 'MDCv6_freqband.dat'
else:
    tableName = 'MDCv6_open_table.dat'
# tablereader returns the following:
#[raInjList, decInjList, fInjList, pulsarNoList, TcohList, sftFileList, PList, asiniList]
tableInfo = tablereader(tableName, observatoryChoice, args)
raInjList = tableInfo[0]
decInjList = tableInfo[1]
fInjList = tableInfo[2]
pulsarNoList = tableInfo[3]
TcohList = tableInfo[4]
sftFileList = tableInfo[5]
PList = tableInfo[6]
asiniList = tableInfo[7]

fSteps = args.fSteps
dfSteps = args.dfSteps

if args.templateSearchClosed:
    print 'Generating Condor files for MDCv6 closed search'
    for m, pulsarNo in enumerate(pulsarNoList):
        shortNo = str(pulsarNo).strip('[').strip(']').strip("'").zfill(3)
        print 'Looking at band, using templateSearchClosed, for pulsar: ' + shortNo
        jobsPerPulsar = 5
        jobInc = m * jobsPerPulsar
        categorizer(TcohList[m][0], raInjList[m][0], decInjList[m][0], fInjList[m][0], observatoryChoice, "pulsar-" + shortNo, sftFileList[m], jobInc, PList[m], asiniList[m], fSteps, dfSteps, args)
else:
    for m, pulsarNo in enumerate(pulsarNoList):
        shortNo = str(pulsarNo).strip('[').strip(']').strip("'").zfill(3)
        print shortNo
        jobsPerPulsar = fSteps * dfSteps
        if args.skyGrid:
            jobsPerPulsar = jobsPerPulsar * args.skyGrid**2
        if args.singleBand:
        # Just one pulsar
            singleBandCode = str(args.singleBand).zfill(3)
            if shortNo == singleBandCode:
                print 'Looking at single band for pulsar: ' + shortNo
                jobInc = 0
                categorizer(TcohList[m][0], raInjList[m][0], decInjList[m][0], fInjList[m][0], observatoryChoice, "pulsar-" + shortNo, sftFileList[m], jobInc, PList[m], asiniList[m], fSteps, dfSteps, args)
        elif args.onlyOne:
        # Again, just one pulsar
            singleCode = str(args.onlyOne).zfill(3)
            if shortNo == singleCode:
                print 'Looking at this pulsar only: ' + shortNo
                jobInc = 0
                categorizer(TcohList[m][0], raInjList[m][0], decInjList[m][0], fInjList[m][0], observatoryChoice, "pulsar-" + shortNo, sftFileList[m], jobInc, PList[m], asiniList[m], fSteps, dfSteps, args)
        elif args.templateSearchOpen:
            print 'Looking at band, using templateSearchOpen, for pulsar: ' + shortNo
            jobsPerPulsar = 5
            jobInc = m * jobsPerPulsar
            categorizer(TcohList[m][0], raInjList[m][0], decInjList[m][0], fInjList[m][0], observatoryChoice, "pulsar-" + shortNo, sftFileList[m], jobInc, PList[m], asiniList[m], fSteps, dfSteps, args)
        elif args.templateSearch:
            if args.templateSearchOpen == false:
            # Using the templateSearch option to scan a 5 Hz band
                templateSearchCode = str(args.templateSearch).zfill(3)
                if shortNo == templateSearchCode:
                    print 'Looking at single band, using templateSearch, for pulsar: ' + shortNo
                    jobInc = 0
                    categorizer(TcohList[m][0], raInjList[m][0], decInjList[m][0], fInjList[m][0], observatoryChoice, "pulsar-" + shortNo, sftFileList[m], jobInc, PList[m], asiniList[m], fSteps, dfSteps, args)
        else:
        # All the pulsars
            if args.noiseTest:
                jobsPerPulsar = int(float(5) / args.noiseTest + 1) * dfSteps
            jobInc = m * jobsPerPulsar
            categorizer(TcohList[m][0], raInjList[m][0], decInjList[m][0], fInjList[m][0], observatoryChoice, "pulsar-" + shortNo, sftFileList[m], jobInc, PList[m], asiniList[m], fSteps, dfSteps, args)
