#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""

fscanDriver.py - Driver script for calling other code to generates SFTs and turns these into plots showing spectra.



"""

__author__ = 'Rejean Dupuis <rejean@caltech.edu> & Greg Mendell<gmendell@ligo-wa.caltech.edu> & Colin Gill <c.gill@astro.gla.ac.uk>'
__date__ = '$Date$'
__version__ = '$Revision$'

# REVISIONS:
# 03/02/2009 gam; implement -g --segment-file option; if not given run on Science,Injections times, if -g ALL given then use start and end times as segment.
# 03/02/2009 gam; if createSFTs then splice in the SFT DAG.
# 06/29/2009 gam; Add option to include matlab path, and run_plotSpecAvgOutput.sh rather than plotSpecAvgOutput.
# 06/29/2009 gam; Include links to _2 plots, when making a comparisons with other fscans.
# 06/29/2009 gam; Fix options and printing of extra debugging info.
# 06/29/2009 gam; use ligolw_segment_query instead of LSCsegFind.
# 08/25/2009 cg; added lines 643 to 652 so that a single frequency band of less than 10 Hz can be run.
# 10/07/2009 gam; Add -I, --intersect-data option to run ligo_data_find with the --show-times option to find times data exist, and use LIGOtools segexpr to intersect this with the segments.
# 10/07/2009 gam; Add -t, --segment-type option to give segment type to use with ligolw_segment_query if segment file is not given (default is IFO:DMT-SCIENCE:1)
# 12/02/2009 gam; Change checked_freqBand to sft_freqBand, and make sure there are some extra bins in the SFTs to avoid problems at endFreq, since lalapps_spec_avg works on [startFreq,endFreq], not [startFreq,endFreq).
#25/02/2010 cg; Added extra options to pass the path of ephemerides for source frequency calculations done in spec_avg.c
#26/04/2010 cg; replaced link to pdf with link to .png files.

# import standard modules and append the lalapps prefix to the python path
import sys, os
import getopt, re, string
import tempfile
import math #cg, so I can use the ceiling functionS
#import ConfigParser
#sys.path.append('')

# import the modules we need to build the pipeline
#from glue import pipeline
#import strain

#
# USAGE FUNCTION
#
def usage():
  msg = """\

Driver script for calling other code to generates SFTs and turns these into plots showing spectra.
  
Usage: [options]

  -h, --help                 display this message
  
  -R, --run                  For a trial run do NOT give this option!
                             When given this code will run condor_submit_dag!
                             Otherwise this script generates the .dag file and then stops!

  -s, --analysis-start-time  GPS start time of data from which to generate SFTs and start analysis
  -L, --duration             length (duration) of data to analyze in seconds
  -G, --tag-string           tag string used in names of various files unique to jobs that will run under the DAG
  -d, --input-data-type      input data type for use with the LSCdataFind --type option
  -x, --extra-datafind-time  (optional) extra time to +/- from/to start/end time used with LSCdataFind (default 256 sec.)
  -M, --datafind-match       (optional) string to use with the LSCdataFind --match option
  -k, --filter-knee-freq     (optional) high pass filter knee frequency used on time domain data before generating SFTs (default = 40 Hz)
  -T, --time-baseline        time baseline of SFTs  (e.g., 60 or 1800 seconds)
  -p, --sft-path             path to SFTs (either already existing or where to output them)
  -c, --freq-res	     the frequency resolution of the spectragrams produced, ***N.B.*** this must be a multiple of 1/T.
  -v, --sft-version          (optional) SFT version number (1 or 2; default is 1)
  -C  --create-sfts          (optional) create the SFTs !!! (/tmp will be appended to the sft-path and SFTs will be generated there!)
  -g, --segment-file         (optional) used only if creating sfts; gives alternative file with segments to use, otherwise ligolw_segment_query is used to find Science segments; if -g ALL is given then [start_time, start_time + duration) is used.
  -t, --segment-type         (optional) segment type to use with ligolw_segment_query if segment file is not given (default is IFO:DMT-SCIENCE:4).
  -o, --sub-log-path         (optional) path to log files given in .sub files (default is $PWD/logs; this directory must exist and usually should be under a local file system.)
  -N, --channel-name         name of input time-domain channel to read from frames
  -i, --ifo                  (optional) ifo to use with ligolw_segment_query and MakeSFTDAG: e.g., H1, H2, L1, G1; PEM channels can start with H0, L0, or G0 (default: use start of channel name)
  -I, --intersect-data       (optional) Run ligo_data_find with the --show-times option to find times data exist, and use LIGOtools segexpr to intersect this with the segments.
  -u, --frame-struct-type    (optional) string specifying the input frame structure and data type. Must begin with ADC_ or PROC_ followed by REAL4, REAL8, INT2, INT4, or INT8; default: ADC_REAL4; -H is the same as PROC_REAL8.
  -F, --start-freq           (optional) start frequency of the SFTs (default is 48 Hz).
  -B, --band                 (optional) frequency band of the SFTs (default is 100 Hz).
  -b, --sub-band             (optional) divide frequency band into sub bands of this size (default is 10 Hz).
  -O, --plot-output-path     (optional) if given then Matlab jobs run and put output plots and data in this directory.
  -m, --matlab-path          (optional) Path to matlab installation, e.g., /ldcg/matlab_r2008a. (Required if --plot-output-path is given.)
  -W, --html-filename        (optional) path and filename of output html file that displays output plots and data.
  -r  --html-reference-path  (optional) path to already existing reference output plots and data for display on html page.
  -e  --html-ref-ifo-epoch   (optional) string that give ifo and GPS start and end times of reference plots, e.g.,H1_840412814_845744945.
  -q  --threshold-snr        (optional) if > 0 and a reference is given, then look for coincident lines with the reference spectra above this threshold.
  -y  --coincidence-deltaf   (optional) if a reference and threshold on snr is given, then use this as the coincidence window on frequency.
  -X  --misc-desc            (optional) misc. part of the SFT description field in the filename (also used if -D option is > 0).
  -H, --use-hot              (optional) input data is from h(t) calibrated frames (h of t = hot!) (0 or 1).
  -J  --psrInput             (optional) Input the name and path of tempo .par file used to overplot a source's freq onto specgrams
  -j  --psrEphemeris         (optional) Input the name and path of ephemeris file used to overplot a source's freq onto specgrams
  -E  --earthFile            (optional) Input the name and path of the .dat file for the earh.  Only used if source details are given.
  -z  --sunFile              (optional) Input the name and path of the .dat file for the Sun.  Only used if source details are given.
  
"""
  print >> sys.stdout, msg

#letters left for short options, a f, n, A, K, Q, U, V, Y.

################################
# MAIN CODE START HERE 
#

####################################
# PARSE COMMAND LINE OPTIONS 
#
shortop = "s:L:G:d:x:M:k:T:p:f:o:N:i:w:P:u:v:c:F:B:b:O:m:W:r:e:q:y:D:X:g:t:l:J:j:E:z:hSHZCRI" #cg; shortop is a string.
longop = [  #cg; longopt is a list
  "help",
  "analysis-start-time=",
  "duration=",
  "dag-file=",
  "tag-string=",
  "input-data-type=",
  "extra-datafind-time=",
  "datafind-match=",
  "filter-knee-freq=",
  "time-baseline=",
  "sft-path=",
  "freq-res=",
  "create-sfts=",
  "log-path=",
  "sub-log-path=",
  "channel-name=",
  "ifo=",
  "intersect-data",
  "frame-struct-type=",
  "window-type=",
  "overlap-fraction=",
  "sft-version=",
  "comment-field=",
  "start-freq=",
  "band=",
  "sub-band=",
  "plot-output-path=",
  "matlab-path=",
  "html-filename=",
  "html-reference-path=",
  "html-ref-ifo-epoch=",
  "threshold-snr=",
  "coincidence-deltaf=",
  "make-gps-dirs=",
  "misc-desc=",
  "max-num-per-node=",
  "segment-file=",
  "segment-type=",
  "min-seg-length=",
  "use-single=",  
  "use-hot",
  "make-tmp-file",
  "run",
  "psrInput=",
  "psrEphmeris=",
  "earthFile=",
  "sunFile=",
  ]

try:
  opts, args = getopt.getopt(sys.argv[1:], shortop, longop)
except getopt.GetoptError:
  usage()
  sys.exit(1)

#############################################
# INITIALIZE DEFAULT VALUES AND READ IN OPTS
#
analysisStartTime = None
analysisEndTime = None
duration = None
tagString = None
inputDataType = None
extraDatafindTime = 256L
datafindMatch = None
filterKneeFreq = 40
timeBaseline = None
pathToSFTs = None
createSFTs = False
cachePath = "cache"
logPath = "logs"
subLogPath = "logs"
channelName = None
segIFO = None
intersectData = False
makeSFTIFO = None
frameStructType = None
windowType = 1
overlapFraction = 0.0
sftVersion = 1
commentField = None
startFreq = 48.0
freqBand = 50.0
freqSubBand = 10.0
plotOutputPath = None
matlabPath = None
makeMatlabPlots = False
htmlFilename = None
htmlReferenceDir = None
htmlRefIFOEpoch = None
thresholdSNR = -1
coincidenceDeltaF = 0
makeGPSDirs = 0
miscDesc = None
maxNumPerNode = 1
maxLengthAllJobs = None
segmentFile = None
segmentType = None
minSegLength = 0L
useSingle = False
useHoT = False
makeTmpFile = False
runCondorSubmitDag = False
freqRes = 0.2
psrInput = None
psrEphemeris = None
earthFile = None
sunFile = None

for o, a in opts:
  if o in ("-h", "--help"):
    usage()
    sys.exit(0)
  elif o in ("-s", "--analysis-start-time"):
    analysisStartTime = long(a)
  elif o in ("-L", "--duration"):
    duration = long(a)
  elif o in ("-G", "--tag-string"):
    tagString = a
  elif o in ("-d", "--input-data-type"):
    inputDataType = a
  elif o in ("-x", "--extra-datafind-time"):
    extraDatafindTime = long(a)
  elif o in ("-M", "--datafind-match"):
    datafindMatch = a
  elif o in ("-k", "--filter-knee-freq"):
    filterKneeFreq = int(a)
  elif o in ("-T", "--time-baseline"):
    timeBaseline = long(a)
  elif o in ("-p", "--sft-path"):
    pathToSFTs = a
  elif o in ("-c", "--freq-res"):
    freqRes = float(a)
  elif o in ("-C", "--create-sfts"):
    createSFTs = True
  elif o in ("-o", "--sub-log-path"):
    subLogPath = a
  elif o in ("-N", "--channel-name"):
    channelName = a
  elif o in ("-i", "--ifo"):
    segIFO = a
  elif o in ("-I", "--intersect-data"):
    intersectData = True
  elif o in ("-u", "--frame-struct-type"):
    frameStructType = a
  elif o in ("-w", "--window-type"):
    windowType = int(a)
  elif o in ("-P", "--overlap-fraction"):
    overlapFraction = float(a)
  elif o in ("-v", "--sft-version"):
    sftVersion = int(a)
  elif o in ("-F", "--start-freq"):
    startFreq = long(a)
  elif o in ("-B", "--band"):
    freqBand = long(a)
  elif o in ("-b", "--sub-band"):
    freqSubBand = long(a)
  elif o in ("-O", "--plot-output-path"):
    plotOutputPath = a
  elif o in ("-m", "--matlab-path"):
    matlabPath = a
  elif o in ("-W", "--html-filename"):
    htmlFilename = a
  elif o in ("-r", "--html-reference-path"):
    htmlReferenceDir = a
  elif o in ("-e", "--html-ref-ifo-epoch"):
    htmlRefIFOEpoch = a
  elif o in ("-q", "--threshold-snr"):
    thresholdSNR = long(a)
  elif o in ("-y", "--coincidence-deltaf"):
    coincidenceDeltaF = long(a)
  elif o in ("-D", "--make-gps-dirs"):
    makeGPSDirs = int(a)
  elif o in ("-X", "--misc-desc"):
    miscDesc = a
  elif o in ("-m", "--max-num-per-node"):
    maxNumPerNode = long(a)
  elif o in ("-L", "--max-length-all-jobs"):
    maxLengthAllJobs = long(a)
  elif o in ("-g", "--segment-file"):
    segmentFile = a
  elif o in ("-t", "--segment-type"):
    segmentType = a
  elif o in ("-l", "--min-seg-length"):
    minSegLength = long(a)
  elif o in ("-S", "--use-single"):
    useSingle = True    
  elif o in ("-H", "--use-hot"):
    useHoT = True
  elif o in ("-R", "--run"):
    runCondorSubmitDag = True
  elif o in ("-Z", "--make-tmp-file"):
    makeTmpFile = True
  elif o in ("-J", "--psrInput"):
    psrInput = a
  elif o in ("-j", "--psrEphmeris"):
    psrEphemeris = a
  elif o in ("-E", "--earthFile"):
    earthFile = a
  elif o in ("-z", "--sunFile"):
    sunFile = a
  else:
    print >> sys.stderr, "Unknown option:", o
    usage()
    sys.exit(1)

#############################################
# VET OPTIONS
#    
if not analysisStartTime:
  print >> sys.stderr, "No analysisStartTime specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not duration:
  print >> sys.stderr, "No duration specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)
 
analysisEndTime = analysisStartTime + duration

if not tagString:
  print >> sys.stderr, "No tag string specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not inputDataType:
  print >> sys.stderr, "No input data type specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if extraDatafindTime < 0L:
  print >> sys.stderr, "Invalid extra datafind time specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)
  
if filterKneeFreq < 0:
  print >> sys.stderr, "No filter knee frequency specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not timeBaseline:
  print >> sys.stderr, "No time baseline specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not pathToSFTs:
  print >> sys.stderr, "No output SFT path specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)
  
if not cachePath:
  print >> sys.stderr, "No cache path specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not logPath:
  print >> sys.stderr, "No log path specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not subLogPath:
  print >> sys.stderr, "No sub log path specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not channelName:
  print >> sys.stderr, "No channel name specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (windowType != 0) and (windowType != 1) and (windowType != 2) and (windowType != 3):
  print >> sys.stderr, "Invalid window type specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (overlapFraction < 0.0) or (overlapFraction >= 1.0):
  print >> sys.stderr, "Invalid make overlap fraction specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (sftVersion != 1) and (sftVersion != 2):
  print >> sys.stderr, "Invalid SFT version specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (startFreq < 0.0):
  print >> sys.stderr, "Invalid start freq specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (freqBand < 0.0):
  print >> sys.stderr, "Invalid band specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (freqSubBand < 0.0):
  print >> sys.stderr, "Invalid sub-band specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (makeGPSDirs < 0) or (makeGPSDirs > 10):
  print >> sys.stderr, "Invalid make gps dirs specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (freqRes < 0.001) or (freqRes > 100.0):
  print >> sys.stderr, "Specify a frequency resoltuion between 100 and 0.001"
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if not maxNumPerNode:
  print >> sys.stderr, "No maximum number of SFTs per node specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (plotOutputPath == None):
   makeMatlabPlots = False
else:
   makeMatlabPlots = True
   if (matlabPath == None):
     print >> sys.stderr, "No matlab path specified."
     print >> sys.stderr, "Use --help for usage details."
     sys.exit(1)

# try and make a directory to store the cache files and job logs
try: os.makedirs(logPath)
except: pass
try: os.makedirs(cachePath)
except: pass
  
# Get site and ifo from channel name:
site = channelName[0]
ifo = channelName[0] + channelName[1]
if not segIFO:
  segIFO = ifo
else:
  # 01/14/07 gam; also send this to MakeSFTDAG using -i option.
  ifo = segIFO
  makeSFTIFO = segIFO

print >> sys.stdout,'\nTHE FSCAN DRIVER SCRIPT HAS STARTED!\n'
  
###################################################
# CHECK IF SFTS NEED TO BE GENERATED
#    
if (createSFTs):
  
  # For safety, add /tmp to the path to avoid overwriting existing SFTs.
  pathToSFTs = pathToSFTs + '/tmp'
  print >> sys.stdout,'Will generate SFTs in %s \n' % pathToSFTs
  try: os.makedirs(pathToSFTs)
  except: pass

  ###################################################
  # FIND SCIENCE MODE SEGMENTS; RUN ligolw_segment_query
  #
  # 03/02/2009 gam; in default case change Science to Science,Injection; add segment file options as elif and else options:
  # 06/29/2009 gam; use ligolw_segment_query instead of LSCsegFind.
  if not segmentFile:
    segmentFile = 'tmpSegs%stmp.txt' % tagString
    if not segmentType:
       # set default segment type
       segmentType = "%s:DMT-SCIENCE:4" % segIFO
    #segCommand = 'LSCsegFind --type Science,Injection --interferometer %s --gps-start-time %d --gps-end-time %d > %s' % (segIFO,analysisStartTime, analysisEndTime,segmentFile)
    segCommand = "ligolw_segment_query --database --query-segments --include-segments %s --gps-start-time %d --gps-end-time %d | grep -v \"0, 0\" | ligolw_print -t segment:table -c start_time -c end_time -d \" \" > %s" % (segmentType,analysisStartTime,analysisEndTime,segmentFile)
    print >> sys.stdout,"Trying: ",segCommand,"\n"
    try:
      segFindExit = os.system(segCommand)
      if (segFindExit > 0):
         print >> sys.stderr, 'ligolw_segment_query failed: %s \n' % segFindExit
         sys.exit(1)
      else:
         print >> sys.stderr, 'ligolw_segment_query succeeded! \n'
    except:
      print >> sys.stderr, 'ligolw_segment_query failed: %s \n' % segFindExit
      sys.exit(1)
  elif segmentFile == "ALL":
    segmentFile = 'tmpSegs%stmp.txt' % tagString
    segCommand = '/bin/echo %d %d > %s' % (analysisStartTime,analysisEndTime,segmentFile)
    print >> sys.stdout,"Trying: ",segCommand,"\n"
    try:
      segFindExit = os.system(segCommand)
      if (segFindExit > 0):
         print >> sys.stderr, 'Failed: %s \n' % segFindExit
         sys.exit(1)
      else:
         print >> sys.stderr, 'Succeeded! \n'
    except:
      print >> sys.stderr, 'Failed: %s \n' % segFindExit
      sys.exit(1)  
  else:
    # Just continue with the name given on the command line.
    print >> sys.stderr, 'Using segmentFile == %s \n' % segmentFile

  ###################################################
  # Intersect with existing data if intersectData = True
  #

  if intersectData:
    # Get the segments the data exist from ligo_data_find
    dataFindSegmentFile = 'tmpDataFindSegs%stmp.txt' % tagString 
    dataFindCommand = "ligo_data_find -s %d -e %d -o %s -t %s -u file --lal-cache --show-times > %s" % (analysisStartTime,analysisEndTime,site,inputDataType,dataFindSegmentFile)
    print >> sys.stdout,"Trying: ",dataFindCommand,"\n"
    try:
      dataFindExit = os.system(dataFindCommand)
      if (dataFindExit > 0):
         print >> sys.stderr, 'ligo_data_find failed: %s \n' % dataFindExit
         sys.exit(1)
      else:
         print >> sys.stderr, 'ligo_data_find succeeded! \n'
    except:
      print >> sys.stderr, 'ligo_data_find failed: %s \n' % dataFindExit
      sys.exit(1)

    # Intersect the segments to run on with the segments data exist
    intersectedSegmentFile = 'tmpIntersectedSegs%stmp.txt' % tagString
    segexprCommand = "segexpr \"intersection(%s,%s)\" %s" % (segmentFile,dataFindSegmentFile,intersectedSegmentFile)
    print >> sys.stdout,"Trying: ",segexprCommand,"\n"
    try:
      segexprExit = os.system(segexprCommand)
      if (segexprExit > 0):
         print >> sys.stderr, 'segexpr failed: %s \n' % segexprExit
         sys.exit(1)
      else:
         print >> sys.stderr, 'segexpr succeeded! \n'
    except:
      print >> sys.stderr, 'segexpr failed: %s \n' % segexprExit
      sys.exit(1)

    # set the segmentFile equal to the intersection of the segments to run on with the segments data exist
    segmentFile = intersectedSegmentFile
  # End if intersectData
    
  ###################################################
  # CHECK THE SEGMENT FILE
  #
  segList = [];
  minSegLength = timeBaseline
  adjustSegExtraTime = True
  try:
    for line in open(segmentFile):
        try: 
            splitLine = line.split();
            try: 
                oneSeg = [];
                oneSeg.append(long(splitLine[0]));
                oneSeg.append(long(splitLine[1]));
                if ((oneSeg[1] - oneSeg[0]) >= minSegLength):
                    segList.append(oneSeg)
                else:
                    pass
            except:
                pass
        except:
            pass
    # End for line in open(segmentFile)
    if (len(segList) < 1):
       print >> sys.stderr, "No segments found in segment file: %s. \n" % segmentFile
       sys.exit(1)
  except:
    print >> sys.stderr, "Error reading or parsing segment file: %s. \n" % segmentFile
    sys.exit(1)
  
  ###################################################
  # MAKE THE .dag FILE; RUN MakeSFTDAG
  #
  sftDAGFile = 'tmpSFTDAG%stmp.dag' % tagString
  # Make sure there are some extra bins in the SFTs to avoid problems at endFreq, since lalapps_spec_avg works on [startFreq,endFreq], not [startFreq,endFreq). 
  if (freqBand < 2):
     sft_freqBand = 2
  else:
     sft_freqBand = freqBand + 1;
  makeDAGCommand = 'MakeSFTDAG -f %s -G %s -d %s -x %d -k %d -T %d -F %d -B %d -p %s -N %s -m 1 -o %s -X %s -Z -g %s -v %d' % (sftDAGFile,tagString,inputDataType,extraDatafindTime,filterKneeFreq,timeBaseline,startFreq,sft_freqBand,pathToSFTs,channelName,subLogPath,miscDesc,segmentFile,sftVersion)
  if (useHoT):
     makeDAGCommand = makeDAGCommand + ' -H'
  if (makeSFTIFO != None):
     makeDAGCommand = makeDAGCommand + ' -i %s' % makeSFTIFO
  if (frameStructType != None):
     makeDAGCommand = makeDAGCommand + ' -u %s' % frameStructType
  print >> sys.stdout,"Trying: ",makeDAGCommand,"\n"
  try:
    makeDAGExit = os.system(makeDAGCommand)
    if (makeDAGExit > 0):
       print >> sys.stderr, 'MakeSFTDAG failed: %s \n' % makeDAGExit
       sys.exit(1)
    else:
       print >> sys.stderr, 'MakeSFTDAG succeeded! \n'
  except:
    print >> sys.stderr, 'MakeSFTDAG failed: %s \n' % makeDAGExit
    sys.exit(1)
else:
  # else if not createSFTs the SFTs already exist, so just continue
  sftDAGFile = None
  print >> sys.stdout,'Will use SFTs already in %s \n' % pathToSFTs


#####################################################
# WRITE A SUPER DAG WHICH RUNS EVERYTHING
#
dagFileName = 'tmpSUPERDAG%stmp.dag' % tagString
dagFileName = 'SUPERDAG%stmp.dag' % tagString
dagFID = file(dagFileName,'w')

# 03/02/2009 gam; use SPLICE to insert SFT DAG
#sftDAGSUBJobName = None
spliceSFTDAGName = None
if (createSFTs):
  # IF GENERATING SFTS ADD THIS TO THE SUPER DAG WHICH RUNS EVERYTHING
  # First write a submit file for summitting a dag to condor
  # 03/02/2009 gam; comment out next lines; no longer need to create condor sub file:
  #condorDAGSUBFileFID = file('condorDAGSUBFile.sub','w')
  #condorDAGSUBFileLogFile = subLogPath + '/' + 'condorDAGSUBFile_' + dagFileName + '.log'
  #condorDAGSUBFileFID.write('universe = scheduler\n')
  #condorDAGSUBFileFID.write('executable = $ENV(CONDOR_LOCATION)/bin/condor_dagman\n')
  #condorDAGSUBFileFID.write('getenv = True\n')
  #condorDAGSUBFileFID.write('arguments = $(argList)\n')
  #condorDAGSUBFileFID.write('log = %s\n' % condorDAGSUBFileLogFile)
  #condorDAGSUBFileFID.write('error = %s/condorDAGSUBFile_$(tagstring).err\n' % logPath)
  #condorDAGSUBFileFID.write('output = %s/condorDAGSUBFile_$(tagstring).out\n' % logPath)
  #condorDAGSUBFileFID.write('notification = never\n')
  #condorDAGSUBFileFID.write('queue 1\n')
  #condorDAGSUBFileFID.close
  ## Now add a jobs to run the sft dag to the super dag
  #sftNodeCount = 0L
  #sftDAGSUBJobName = 'runSFTDAGSUB_%i' % sftNodeCount
  #dagFID.write('JOB %s condorDAGSUBFile.sub\n' % sftDAGSUBJobName)
  #dagFID.write('RETRY %s 1\n' % sftDAGSUBJobName)
  #sftDAGSUBFileLogFile = subLogPath + '/' + 'sftDAGSUBFile_' + dagFileName + '.log'
  #argList = '-f -l . -Debug 3 -Lockfile %s.lock -Condorlog %s -Dag %s -Rescue %s.rescue -MaxJobs 40 -MaxPre 1 -MaxPost 1' % (sftDAGFile, sftDAGSUBFileLogFile, sftDAGFile, sftDAGFile)
  #tagStringOut = '%s_%i' % (tagString, sftNodeCount)
  #dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(sftDAGSUBJobName,argList,tagStringOut))
  
  # 03/02/2009 gam; instead of above, add splice in SFT DAG:
  sftNodeCount = 0L
  spliceSFTDAGName = 'spliceSFTDAG_%i' % sftNodeCount
  dagFID.write('SPLICE %s %s\n' % (spliceSFTDAGName, sftDAGFile))

# CREATE A CONDOR .sub FILE TO RUN lalapps_spec_avg
spectrumAverageFID = file('spectrumAverage.sub','w') #cg; open then file for writing
spectrumAverageLogFile = subLogPath + '/' + 'spectrumAverage_' + dagFileName + '.log'
spectrumAverageFID.write('universe = vanilla\n') #cg; write these followig lines into the file.
spectrumAverageFID.write('executable = $ENV(SPECAVG_PATH)/lalapps_spec_avg\n')
spectrumAverageFID.write('arguments = $(argList)\n')
spectrumAverageFID.write('log = %s\n' % spectrumAverageLogFile)
spectrumAverageFID.write('error = %s/spectrumAverage_$(tagstring).err\n' % logPath)
spectrumAverageFID.write('output = %s/spectrumAverage_$(tagstring).out\n' % logPath)
spectrumAverageFID.write('notification = never\n')
spectrumAverageFID.write('queue 1\n')
spectrumAverageFID.close()

# MAKE A SUBMIT FILE FOR RUNNING THE MATLAB DRIVER SCRIPT
if (makeMatlabPlots):
  runMatlabScriptFID = file('runMatlabPlotScript.sub','w')
  runMatlabScriptLogFile = subLogPath + '/' + 'runMatlabPlotScript_' + dagFileName + '.log'
  runMatlabScriptFID.write('universe = local\n')
  # Run compiled version plotSpecAvgOutput.m:
  # runMatlabScriptFID.write('executable = $ENV(PLOTSPECAVGOUTPUT_PATH)/plotSpecAvgOutput\n'); # 06/29/09 gam; changed to run_plotSpecAvgOutput.sh. 
  runMatlabScriptFID.write('executable = $ENV(PLOTSPECAVGOUTPUT_PATH)/run_plotSpecAvgOutput.sh\n')
  runMatlabScriptFID.write('getenv = True\n')
  runMatlabScriptFID.write('arguments = $(argList)\n')
  runMatlabScriptFID.write('log = %s\n' % runMatlabScriptLogFile)
  runMatlabScriptFID.write('error = %s/runMatlabPlotScript_$(tagstring).err\n' % logPath)
  runMatlabScriptFID.write('output = %s/runMatlabPlotScript_$(tagstring).out\n' % logPath)
  runMatlabScriptFID.write('notification = never\n')
  runMatlabScriptFID.write('queue 1\n')
  runMatlabScriptFID.close()

#cg; Creates the html file, or at least the initial parts of it.
if (htmlFilename != None):
  htmlFID = file(htmlFilename,'w')#opens the html file specified in the command line for writing, is having trouble for some reason.
  htmlFID.write('<html>\n')
  htmlFID.write('<head>\n')
  htmlFID.write('<meta content="text/html; charset=ISO-8859-1" http-equiv="content-type">\n')
  htmlFID.write('<title>FSCANPLOTS</title>\n')
  htmlFID.write('</head>\n')
  htmlFID.write('<body>\n')
  htmlFID.write('<div style="text-align: center;">\n')
  htmlFID.write('<h1>FSCAN PLOTS</h1>\n')
  htmlFID.write('</div>\n')
  htmlFID.write('<br>\n')
  if (htmlReferenceDir != None) and (thresholdSNR > 0):#If we want to compare to a reference, then lines below will create an extra page, and link to it...
    htmlLinesFilenameList = os.path.split(htmlFilename)
    htmlLinesFilenameShort = 'Lines_%s' % htmlLinesFilenameList[1]
    htmlLinesFilename = os.path.join(htmlLinesFilenameList[0],htmlLinesFilenameShort)
    htmlLinesFID = file(htmlLinesFilename,'w')#create a file whose name is stored in htmlLinesFilename, open it for writing.
    htmlLinesFID.write('<html>\n')#write html into file, and go to next line.
    htmlLinesFID.write('<head>\n')
    htmlLinesFID.write('<meta content="text/html; charset=ISO-8859-1" http-equiv="content-type">\n')
    htmlLinesFID.write('<title>FSCANCOINCIDENTLINES</title>\n')
    htmlLinesFID.write('</head>\n')
    htmlLinesFID.write('<body>\n')
    htmlLinesFID.write('<div style="text-align: center;">\n')
    htmlLinesFID.write('<h1>FSCAN COINCIDENT LINES WITH SNR >= %d</h1>\n' % thresholdSNR)
    htmlLinesFID.write('</div>\n')
    htmlLinesFID.write('<br>\n')  
    # Add a link from the main page to this file:
    htmlFID.write('<div style="text-align: center;">\n')    
    htmlFID.write('Click <a href="%s">here</a> to get a list of coincident lines (or click on the link to "Coincident Lines" below each plot on the left).<br>\n' % htmlLinesFilenameShort)
    htmlFID.write('</div>\n')
    htmlFID.write('<br>\n')
  htmlFID.write('<div style="text-align: center;">\n')
  htmlFID.write('<h3>Click on a plot to see the .png source file. Click on a link below a plot to get the SFT timestamps, the spectrogram data file, frequency vs power  and SNR data file.</h3>\n')
  htmlFID.write('</div>\n')
  htmlFID.write('<br>\n')  
  htmlFID.write('<form>\n')
  htmlFID.write('Old <input type = "radio" name = "rad1" onclick = "showDiv()" checked>\n')
  htmlFID.write('New <input type = "radio" name = "rad1" onclick = "showDiv()">\n')
  htmlFID.write('<br><br>\n')
  htmlFID.write('<div id = "div1" style="display:none">\n')
  htmlFID.write('<table style="width: 100%; text-align: left;" border="1" cellpadding="2" cellspacing="2">\n')
  htmlFID.write('<tbody>\n')


if (freqSubBand  > freqBand):
  freqSubBand = freqBand

endFreq = startFreq + freqBand
thisStartFreq = startFreq
thisEndFreq = thisStartFreq + freqSubBand
nodeCount = 0L

# Next lines should not be needed because sft_freqBand goes beyond endFreq above
#if (thisEndFreq == endFreq):
#  endFreq=endFreq+1

# Note that sft_freqBand goes beyond endFreq above, so problem that lalapps_spec_avg
# works on [startFreq,endFreq], not [startFreq,endFreq), should now be avoided.
# Thus, can use "while (thisEndFreq <= endFreq)" rather than "while (thisEndFreq < endFreq)".
while (thisEndFreq <= endFreq):

  #print >> sys.stdout,"startfreq: ",startFreq,"\n"
  #print >> sys.stdout,"endfreq: ",endFreq,"\n"
  #print >> sys.stdout,"thisStartfreq: ",thisStartFreq,"\n"
  #print >> sys.stdout,"thisEndfreq: ",thisEndFreq,"\n"
  
  #check that the freq resolution does not exceed that of the raw sft freq resolution, if it does replace by max resolution. Also need to check that freqres is an integer multiple of the raw sft freqres, if not replace it with one that is.
  
  NumBinsAvg = math.floor(freqRes*timeBaseline);
  freqRes=(1/float(timeBaseline))*NumBinsAvg#sets the freqRes so its an multiple of the raw sft freqres.
  if (freqRes < (1/float(timeBaseline))):#makes sure that this is not less than raw sft freqres, it mght be zero.
      freqRes = (1/float(timeBaseline))#if it is zero, just make it so that its the same res as the raw sft res
      print >> sys.stderr, 'freqRes changed to: %f \n' % freqRes#let user know its been changed,
  effTBaseFull = timeBaseline
  effTBase = 1/freqRes
  #end of freqres checking

  fRange = thisEndFreq - thisStartFreq

  #This section of code writes the arguments for spec_avg.c
  #------------------------------------------------------------
  specAvgJobName = 'SpecAvg_%i' % nodeCount
  dagFID.write('JOB %s spectrumAverage.sub\n' % specAvgJobName)
  dagFID.write('RETRY %s 0\n' % specAvgJobName)
  #the variable names in the arglist below are the names of the variables in spec_avg.c not in this code
  if ((psrInput != None) and (sftVersion == 2)):
     argList = '--startGPS %d --endGPS %d --IFO %s --fMin %d --fMax %d --freqRes %f  --timeBaseline %d --SFTs %s/*.sft  --psrInput %s --psrEphemeris %s --earthFile %s --sunFile %s' % (analysisStartTime,analysisEndTime,ifo,thisStartFreq,thisEndFreq,freqRes,timeBaseline,pathToSFTs,psrInput, psrEphemeris, earthFile, sunFile)
  elif ((psrInput == None) and (sftVersion == 2)):
     argList = '--startGPS %d --endGPS %d --IFO %s --fMin %d --fMax %d --freqRes %f  --timeBaseline %d --SFTs %s/*.sft' % (analysisStartTime,analysisEndTime,ifo,thisStartFreq,thisEndFreq,freqRes,timeBaseline,pathToSFTs)
  elif ((psrInput != None) and (sftVersion != 2)):
     argList = '--startGPS %d --endGPS %d --IFO %s --fMin %d --fMax %d --freqRes %f --timeBaseline %d --SFTs %s/SFT*   --psrInput %s --psrEphemeris %s --earthFile %s --sunFile %s' % (analysisStartTime,analysisEndTime,ifo,thisStartFreq,thisEndFreq,freqRes,timeBaseline,pathToSFTs,psrInput, psrEphemeris, earthFile, sunFile)
  else:
     argList = '--startGPS %d --endGPS %d --IFO %s --fMin %d --fMax %d --freqRes %f --timeBaseline %d --SFTs %s/SFT*' % (analysisStartTime,analysisEndTime,ifo,thisStartFreq,thisEndFreq,freqRes,timeBaseline,pathToSFTs)
  tagStringOut = '%s_%i' % (tagString, nodeCount)  
  dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(specAvgJobName,argList,tagStringOut))
  if (createSFTs):
    # 03/02/2009 gam; use SPLICE to insert SFT DAG
    #dagFID.write('PARENT %s CHILD %s\n'%(sftDAGSUBJobName,specAvgJobName))
    dagFID.write('PARENT %s CHILD %s\n'%(spliceSFTDAGName,specAvgJobName))

#write arguments to be passed to matlab job
  if (makeMatlabPlots):
    runMatlabPlotScriptJobName = 'runMatlabPlotScript_%i' % nodeCount
    dagFID.write('JOB %s runMatlabPlotScript.sub\n' % runMatlabPlotScriptJobName)
    dagFID.write('RETRY %s 0\n' % runMatlabPlotScriptJobName)
    inputFileName = 'spec_%d.00_%d.00_%s_%d_%d' % (thisStartFreq,thisEndFreq,ifo,analysisStartTime,analysisEndTime)
    outputFileName = '%s/%s' % (plotOutputPath, inputFileName)   

    #work out the best frequency interval to place ticks at.
    
    if (fRange >= 200):
	deltaFTicks = 20
    elif (fRange >= 100):
	deltaFTicks = 10
    elif (fRange >= 50):
	deltaFTicks = 5
    elif (fRange >= 5):
	deltaFTicks = 1
    else:
	deltaFTicks = 0.1
    
    taveFlag = 1
    
    #Test to see if we have pulsar input
    if (psrInput != None):
        pulsar = 1
    else:
        pulsar = 0

    if (htmlReferenceDir != None):
	referenceFileName = '%s/spec_%d.00_%d.00_%s' % (htmlReferenceDir,thisStartFreq,thisEndFreq,htmlRefIFOEpoch)
    else:
        referenceFileName = 'none'
    #argList = '%s %s %s %d %d %d %d %d %d %s' % (inputFileName,outputFileName,channelName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coincidenceDeltaF,referenceFileName); # 06/29/09 gam; matlabPath has to be first argument.
    #argList = '%s %s %s %s %d %d %d %d %d %d %s' % (matlabPath,inputFileName,outputFileName,channelName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coincidenceDeltaF,referenceFileName)
    argList = '%s %s %s %s %d %g %d %d %d %d %d %s' % (matlabPath,inputFileName,outputFileName,channelName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coincidenceDeltaF,pulsar,referenceFileName)
    tagStringOut = '%s_%i' % (tagString, nodeCount)  
    dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(runMatlabPlotScriptJobName,argList,tagStringOut))
    dagFID.write('PARENT %s CHILD %s\n'%(specAvgJobName,runMatlabPlotScriptJobName))
  
  #add stuff to the html file to show the plots produced in matlab.
  #-----------------------------------------------------------------
  if (htmlFilename != None):
    inputFileName = 'spec_%d.00_%d.00_%s_%d_%d' % (thisStartFreq,thisEndFreq,ifo,analysisStartTime,analysisEndTime)
    htmlFID.write('  <tr>\n')#tr adds a new table row.
    htmlFID.write('  <td style="vertical-align: top;">\n')
    htmlFID.write('    <a href="%s.png"><img alt="" src="%s.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    <a href="%s_2.png"><img alt="" src="%s_2.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('  </td>\n')#ends the current data cell.
    if (htmlReferenceDir != None):
        referenceFileName = '%s/spec_%d.00_%d.00_%s' % (htmlReferenceDir,thisStartFreq,thisEndFreq,htmlRefIFOEpoch)
        htmlFID.write('  <td style="vertical-align: top;">\n')
        htmlFID.write('    <a href="%s.png"><img alt="" src="%s.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('    <a href="%s_2.png"><img alt="" src="%s_2.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('  </td>\n')     
    htmlFID.write('  </tr>\n')
    # Topten lines code
    #htmlFID.write('  <tr>\n')
    #htmlFID.write('  <td style="vertical-align: top;">\n')
    #htmlFID.write('<object data="%s_topten.txt" type="text/plain" style="width: 620px; height: 450px"></object>' % inputFileName)
    #htmlFID.write('  </td\n')
    #htmlFID.write('  </tr>\n')
    #end of topten lines.
    htmlFID.write('  <tr>\n')#open new row in table for lines files.
    htmlFID.write('  <td style="vertical-align: top;">\n')
    htmlFID.write('    SFT Timestamps: <a href="%s_timestamps">%s_timestamps</a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    Spectrogram data: <a href="%s">%s</a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    Freq. vs Power, SNR: <a href="%s.txt">%s.txt</a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    Freq. vs Power, SNR (Sorted): <a href="%s_sorted.txt">%s_sorted.txt</a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    List of found combs : <a href="%s_combs.txt">%s_combs.txt</a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    Kurtosis test output: <a href="%s_kurtosis">%s_kurtosis</a><br>\n' % (inputFileName,inputFileName))
    if (htmlReferenceDir != None) and (thresholdSNR > 0):
        #--------------------------
	#coincident lines
	htmlFID.write('    Coincident Lines: <a href="%s_coincident_lines.txt">%s_coincident_lines.txt</a><br>\n' % (inputFileName,inputFileName))
        htmlLinesFID.write('<br>\n')
        htmlLinesFID.write('<object data="%s_coincident_lines.txt" type="text/plain" style="width: 620px; height: 450px"></object>' % inputFileName)
        htmlLinesFID.write('<br>\n')
	#--------------------------
	#new lines
	htmlFID.write('    New Lines: <a href="%s_new_lines.txt">%s_new_lines.txt</a><br>\n' % (inputFileName,inputFileName))
        htmlLinesFID.write('<br>\n')
        htmlLinesFID.write('<object data="%s_new_lines.txt" type="text/plain" style="width: 620px; height: 450px"></object>' % inputFileName)
        htmlLinesFID.write('<br>\n')
	#---------------------------
	#old lines
	htmlFID.write('    Old Lines: <a href="%s_old_lines.txt">%s_old_lines.txt</a><br>\n' % (inputFileName,inputFileName))
        htmlLinesFID.write('<br>\n')
        htmlLinesFID.write('<object data="%s_old_lines.txt" type="text/plain" style="width: 620px; height: 450px"></object>' % inputFileName)
        htmlLinesFID.write('<br>\n')
	#---------------------------
    htmlFID.write('  </td>\n')
    if (htmlReferenceDir != None):
        htmlFID.write('  <td style="vertical-align: top;">\n')
        htmlFID.write('    SFT Timestamps: <a href="%s_timestamps">%s_timestamps</a><br>\n' % (referenceFileName,referenceFileName))
	htmlFID.write('    Spectrogram data: <a href="%s">%s</a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('    Freq. vs Power, SNR: <a href="%s.txt">%s.txt</a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('    Freq. vs Power, SNR (Sorted): <a href="%s_sorted.txt">%s_sorted.txt</a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('    List of found combs : <a href="%s_combs.txt">%s_combs.txt</a><br>\n' % (inputFileName,inputFileName))
        htmlFID.write('  </td>\n')
    htmlFID.write('  </tr>\n')
  thisStartFreq = thisStartFreq + freqSubBand
  thisEndFreq = thisStartFreq + freqSubBand
  nodeCount = nodeCount + 1

# Close the DAG file
dagFID.close()

if (htmlFilename != None):
  htmlFID.write('</tbody>\n')
  htmlFID.write('</table>\n')
  htmlFID.write('</div>\n')
  htmlFID.write('<div id = "div2" style="display:none">\n')
  htmlFID.write('<table>\n')
  rows =  freqBand/freqSubBand
  lastFreq = startFreq +(rows*freqSubBand)
  for i in range(rows):
    thatLastFreq = lastFreq - (freqSubBand*(i+1))
    thisLastFreq = lastFreq - (freqSubBand*i)
    thisSeg = 'spec_%d.00_%d.00_%s_%d_%d' % (thatLastFreq,thisLastFreq,ifo,analysisStartTime,analysisEndTime)
    htmlFID.write('<tr>\n')
    htmlFID.write('  <td>\n')
    htmlFID.write('    <a href="' + thisSeg + '.png"><img alt="" src="' + thisSeg + '.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n')
    htmlFID.write('  </td>\n')
    htmlFID.write('  <td>\n')
    htmlFID.write('    <a href="' + thisSeg + '_2.png"><img alt="" src="' + thisSeg + '_2.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n')
    htmlFID.write('  </td>\n')
    htmlFID.write('  <td>\n')
    htmlFID.write('    SFT Timestamps: <a href="' + thisSeg + '_timestamps">' + thisSeg + '_timestamps</a><br>')
    htmlFID.write('    Spectrogram data: <a href="' + thisSeg + '">' + thisSeg + '</a><br>')
    htmlFID.write('    Freq. vs Power, SNR: <a href="' + thisSeg + '.txt">' + thisSeg + '.txt</a><br>')
    htmlFID.write('    Freq. vs Power, SNR (Sorted): <a href="' + thisSeg + '_sorted.txt">' + thisSeg + '_sorted.txt</a><br>')
    htmlFID.write('    List of found combs : <a href="%s_combs.txt">%s_combs.txt</a><br>\n' % (inputFileName,inputFileName))
    htmlFID.write('    Kurtosis test output: <a href="' + thisSeg + '_kurtosis">' + thisSeg + '_kurtosis</a><br>')
    htmlFID.write('  </td>\n')
    htmlFID.write('</tr>\n')
  htmlFID.write('</table>\n')
  htmlFID.write('</div>\n')
  htmlFID.write('</form>\n')
  htmlFID.write('<script type = "text/javascript">\n')
  htmlFID.write('function showDiv() {\n')
  htmlFID.write('document.getElementById("div1").style.display="none";\n')
  htmlFID.write('document.getElementById("div2").style.display="none";\n')
  htmlFID.write('if (document.forms[0].rad1[0].checked) {\n')
  htmlFID.write('document.getElementById("div1").style.display="block";\n')
  htmlFID.write('}\n')
  htmlFID.write('if (document.forms[0].rad1[1].checked) {\n')
  htmlFID.write('document.getElementById("div2").style.display="block";\n')
  htmlFID.write('}\n')
  htmlFID.write('}\n')
  htmlFID.write('</script>\n')
  htmlFID.write('<span style="text-decoration: underline;"><br>\n')
  htmlFID.write('</span>\n')
  htmlFID.write('</body>\n')
  htmlFID.write('</html>\n')
  htmlFID.close
  if (htmlReferenceDir != None) and (thresholdSNR > 0):
     htmlLinesFID.write('<span style="text-decoration: underline;"><br>\n')
     htmlLinesFID.write('</span>\n')
     htmlLinesFID.write('</body>\n')
     htmlLinesFID.write('</html>\n')
     htmlLinesFID.close()

###################################################
# SUBMIT THE .dag FILE TO CONDOR; RUN condor_submit_dag
#
runDAGCommand = 'condor_submit_dag %s' % dagFileName
print >> sys.stdout,"Trying: ",runDAGCommand,"\n"
if (runCondorSubmitDag):
   try:
       runDAGExit = os.system(runDAGCommand)
       if (runDAGExit > 0):
          print >> sys.stderr, 'condor_submit_dag failed: %s \n' % runDAGExit
          sys.exit(1)
       else:
          print >> sys.stderr, 'condor_submit_dag succeeded! \n'
   except:
       print >> sys.stderr, 'condor_submit_dag failed: %s \n' % runDAGExit
       sys.exit(1)
else:
   print >> sys.stderr, 'TRIAL RUN ONLY!!! Either submit %s by hand or run this script with the -R or --run option! \n' % dagFileName
   sys.exit(1)

###################################################
# CLEAN UP
#        
# rmOut = os.system('/bin/rm -f %s 1>/dev/null 2>/dev/null' % segmentFile    

sys.exit(0)
