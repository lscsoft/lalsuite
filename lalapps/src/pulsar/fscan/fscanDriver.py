#!/usr/bin/env python
"""

fscanDriver.py - Driver script for calling other code to generates SFTs and turns these into plots showing spectra.


$Id$

"""

__author__ = 'Rejean Dupuis <rejean@caltech.edu> & Greg Mendell<gmendell@ligo-wa.caltech.edu>'
__date__ = '$Date$'
__version__ = '$Revision$'[11:-2]

# REVISIONS:

# import standard modules and append the lalapps prefix to the python path
import sys, os
import getopt, re, string
import tempfile
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
  -v, --sft-version          (optional) SFT version number (1 or 2; default is 1)
  -C  --create-sfts          (optional) create the SFTs !!! (/tmp will be appended to the sft-path and SFTs will be generated there!)
  -o, --sub-log-path         (optional) path to log files given in .sub files (default is $PWD/logs; this directory must exist and usually should be under a local file system.)
  -N, --channel-name         name of input time-domain channel to read from frames
  -i, --ifo                  (optional) ifo to use with LSCsegFind and MakeSFTDAG: e.g., H1, H2, L1, G1; PEM channels can start with H0, L0, or G0 (default: use start of channel name)
  -u, --frame-struct-type    (optional) string specifying the input frame structure and data type. Must begin with ADC_ or PROC_ followed by REAL4, REAL8, INT2, INT4, or INT8; default: ADC_REAL4; -H is the same as PROC_REAL8.
  -F, --start-freq           (optional) start frequency of the SFTs (default is 48 Hz).
  -B, --band                 (optional) frequency band of the SFTs (default is 100 Hz).
  -b, --sub-band             (optional) divide frequency band into sub bands of this size (default is 10 Hz).
  -O, --plot-output-path     (optional) if given then Matlab jobs run and put output plots and data in this directory.
  -W, --html-filename        (optional) path and filename of output html file that displays output plots and data.
  -r  --html-reference-path  (optional) path to already existing reference output plots and data for display on html page.
  -e  --html-ref-ifo-epoch   (optional) string that give ifo and GPS start and end times of reference plots, e.g.,H1_840412814_845744945.
  -q  --threshold-snr        (optional) if > 0 and a reference is given, then look for coincident lines with the reference spectra above this threshold.
  -y  --coincidence-deltaf   (optional) if a reference and threshold on snr is given, then use this as the coincidence window on frequency.
  -X  --misc-desc            (optional) misc. part of the SFT description field in the filename (also used if -D option is > 0).
  -H, --use-hot              (optional) input data is from h(t) calibrated frames (h of t = hot!) (0 or 1).
  
"""
  print >> sys.stdout, msg

################################
# MAIN CODE START HERE 
#

####################################
# PARSE COMMAND LINE OPTIONS 
#
shortop = "s:L:G:d:x:M:k:T:p:o:N:i:w:P:u:v:c:F:B:b:O:W:r:e:q:y:D:X:m:g:l:hSHZCR"
longop = [
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
  "create-sfts=",
  "log-path=",
  "sub-log-path=",
  "channel-name=",
  "ifo=",
  "frame-struct-type=",
  "window-type=",
  "overlap-fraction=",
  "sft-version=",
  "comment-field=",
  "start-freq=",
  "band=",
  "sub-band=",
  "plot-output-path=",
  "html-filename=",
  "html-reference-path=",
  "html-ref-ifo-epoch=",
  "threshold-snr=",
  "coincidence-deltaf=",
  "make-gps-dirs=",
  "misc-desc=",
  "max-num-per-node=",
  "segment-file=",
  "min-seg-length=",
  "use-single=",  
  "use-hot",
  "make-tmp-file",
  "run",
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
makeSFTIFO = None
frameStructType = None
windowType = 1
overlapFraction = 0.0
sftVersion = 1
commentField = None
startFreq = 48.0
freqBand = 100.0
freqSubBand = 10.0
plotOutputPath = None
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
minSegLength = 0L
useSingle = False
useHoT = False
makeTmpFile = False
runCondorSubmitDag = False

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
  elif o in ("-C", "--create-sfts"):
    createSFTs = True
  elif o in ("-o", "--sub-log-path"):
    subLogPath = a
  elif o in ("-N", "--channel-name"):
    channelName = a
  elif o in ("-i", "--ifo"):
    segIFO = a
  elif o in ("-u", "--frame-struct-type"):
    frameStructType = a
  elif o in ("-w", "--window-type"):
    windowType = int(a)
  elif o in ("-P", "--overlap-fraction"):
    overlapFraction = float(a)
  elif o in ("-v", "--sft-version"):
    sftVersion = int(a)
  elif o in ("-c", "--comment-field"):
    commentField = a
  elif o in ("-F", "--start-freq"):
    startFreq = long(a)
  elif o in ("-B", "--band"):
    freqBand = long(a)
  elif o in ("-b", "--sub-band"):
    freqSubBand = long(a)
  elif o in ("-O", "--plot-output-path"):
    plotOutputPath = a
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

if not maxNumPerNode:
  print >> sys.stderr, "No maximum number of SFTs per node specified."
  print >> sys.stderr, "Use --help for usage details."
  sys.exit(1)

if (plotOutputPath == None):
   makeMatlabPlots = False
else:
   makeMatlabPlots = True
  
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
  # FIND SCIENCE MODE SEGMENTS; RUN LSCsegFind
  #
  segmentFile = 'tmpSegs%stmp.txt' % tagString
  segCommand = 'LSCsegFind --type Science --interferometer %s --gps-start-time %d --gps-end-time %d > %s' % (segIFO,analysisStartTime, analysisEndTime,segmentFile)
  print >> sys.stdout,"Trying: ",segCommand,"\n"
  try:
    segFindExit = os.system(segCommand)
    if (segFindExit > 0):
       print >> sys.stderr, 'LSCsegFind failed: %s \n' % segFindExit
       sys.exit(1)
    else:
       print >> sys.stderr, 'LSCsegFind succeeded! \n'
  except:
    print >> sys.stderr, 'LSCsegFind failed: %s \n' % segFindExit
    sys.exit(1)
  
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
  makeDAGCommand = 'MakeSFTDAG -f %s -G %s -d %s -x %d -k %d -T %d -F %d -B %d -p %s -N %s -m 1 -o %s -X %s -Z -g %s -v %d' % (sftDAGFile,tagString,inputDataType,extraDatafindTime,filterKneeFreq,timeBaseline,startFreq,freqBand,pathToSFTs,channelName,subLogPath,miscDesc,segmentFile,sftVersion)
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

sftDAGSUBJobName = None
if (createSFTs):
  # IF GENERATING SFTS ADD THIS TO THE SUPER DAG WHICH RUNS EVERYTHING
  # First write a submit file for summitting a dag to condor
  condorDAGSUBFileFID = file('condorDAGSUBFile.sub','w')
  condorDAGSUBFileLogFile = subLogPath + '/' + 'condorDAGSUBFile_' + dagFileName + '.log'
  condorDAGSUBFileFID.write('universe = scheduler\n')
  condorDAGSUBFileFID.write('executable = $ENV(CONDOR_LOCATION)/bin/condor_dagman\n')
  condorDAGSUBFileFID.write('getenv = True\n')
  condorDAGSUBFileFID.write('arguments = $(argList)\n')
  condorDAGSUBFileFID.write('log = %s\n' % condorDAGSUBFileLogFile)
  condorDAGSUBFileFID.write('error = %s/condorDAGSUBFile_$(tagstring).err\n' % logPath)
  condorDAGSUBFileFID.write('output = %s/condorDAGSUBFile_$(tagstring).out\n' % logPath)
  condorDAGSUBFileFID.write('notification = never\n')
  condorDAGSUBFileFID.write('queue 1\n')
  condorDAGSUBFileFID.close
  # Now add a jobs to run the sft dag to the super dag
  sftNodeCount = 0L
  sftDAGSUBJobName = 'runSFTDAGSUB_%i' % sftNodeCount
  dagFID.write('JOB %s condorDAGSUBFile.sub\n' % sftDAGSUBJobName)
  dagFID.write('RETRY %s 1\n' % sftDAGSUBJobName)
  sftDAGSUBFileLogFile = subLogPath + '/' + 'sftDAGSUBFile_' + dagFileName + '.log'
  argList = '-f -l . -Debug 3 -Lockfile %s.lock -Condorlog %s -Dag %s -Rescue %s.rescue -MaxJobs 40 -MaxPre 1 -MaxPost 1' % (sftDAGFile, sftDAGSUBFileLogFile, sftDAGFile, sftDAGFile)
  tagStringOut = '%s_%i' % (tagString, sftNodeCount)
  dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(sftDAGSUBJobName,argList,tagStringOut))

# CREATE A CONDOR .sub FILE TO RUN lalapps_spec_avg
spectrumAverageFID = file('spectrumAverage.sub','w')
spectrumAverageLogFile = subLogPath + '/' + 'spectrumAverage_' + dagFileName + '.log'
spectrumAverageFID.write('universe = vanilla\n')
spectrumAverageFID.write('executable = $ENV(SPECAVG_PATH)/lalapps_spec_avg\n')
spectrumAverageFID.write('arguments = $(argList)\n')
spectrumAverageFID.write('log = %s\n' % spectrumAverageLogFile)
spectrumAverageFID.write('error = %s/spectrumAverage_$(tagstring).err\n' % logPath)
spectrumAverageFID.write('output = %s/spectrumAverage_$(tagstring).out\n' % logPath)
spectrumAverageFID.write('notification = never\n')
spectrumAverageFID.write('queue 1\n')
spectrumAverageFID.close

# MAKE A SUBMIT FILE FOR RUNNING THE MATLAB DRIVER SCRIPT
if (makeMatlabPlots):
  runMatlabScriptFID = file('runMatlabPlotScript.sub','w')
  runMatlabScriptLogFile = subLogPath + '/' + 'runMatlabPlotScript_' + dagFileName + '.log'
  runMatlabScriptFID.write('universe = vanilla\n')
  # Run compiled version plotSpecAvgOutput.m:
  runMatlabScriptFID.write('executable = $ENV(PLOTSPECAVGOUTPUT_PATH)/plotSpecAvgOutput\n')
  runMatlabScriptFID.write('getenv = True\n')
  runMatlabScriptFID.write('arguments = $(argList)\n')
  runMatlabScriptFID.write('log = %s\n' % runMatlabScriptLogFile)
  runMatlabScriptFID.write('error = %s/runMatlabPlotScript_$(tagstring).err\n' % logPath)
  runMatlabScriptFID.write('output = %s/runMatlabPlotScript_$(tagstring).out\n' % logPath)
  runMatlabScriptFID.write('notification = never\n')
  runMatlabScriptFID.write('queue 1\n')
  runMatlabScriptFID.close

if (htmlFilename != None):
  htmlFID = file(htmlFilename,'w')
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
  if (htmlReferenceDir != None) and (thresholdSNR > 0):
    htmlLinesFilenameList = os.path.split(htmlFilename)
    htmlLinesFilenameShort = 'Lines_%s' % htmlLinesFilenameList[1]
    htmlLinesFilename = os.path.join(htmlLinesFilenameList[0],htmlLinesFilenameShort)
    htmlLinesFID = file(htmlLinesFilename,'w')
    htmlLinesFID.write('<html>\n')
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
  htmlFID.write('<h3>Click on a plot to get a larger pdf version. Click on a link below a plot to get the corresponding timestamp or frequency vs power data file.</h3>\n')
  htmlFID.write('</div>\n')
  htmlFID.write('<br>\n')  
  htmlFID.write('<table style="width: 100%; text-align: left;" border="1" cellpadding="2" cellspacing="2">\n')
  htmlFID.write('<tbody>\n')

# Write lalapps_spec_avg jobs to SUPER DAG:
endFreq = startFreq + freqBand
thisStartFreq = startFreq
thisEndFreq = thisStartFreq + freqSubBand
nodeCount = 0L
while (thisEndFreq < endFreq):
  # Because SFTs are in the band [startFreq,endFreq) but lalapps_spec_avg works on [startFreq,endFreq]
  # we need to avoid the last subBand requested. The user thus should request 1 Hz more than
  # wanted and the last one Hz band will not be used.
  #if (thisEndFreq >= endFreq):
     # Fix off-by-one bin end frequency error; for simplicity remove final 1 Hz
     #thisEndFreq = endFreq - 1
  specAvgJobName = 'SpecAvg_%i' % nodeCount
  dagFID.write('JOB %s spectrumAverage.sub\n' % specAvgJobName)
  dagFID.write('RETRY %s 10\n' % specAvgJobName)
  if (sftVersion == 2):
     argList = '--startGPS %d --endGPS %d --IFO %s --fMin %d --fMax %d --SFTs %s/*.sft' % (analysisStartTime,analysisEndTime,ifo,thisStartFreq,thisEndFreq,pathToSFTs)
  else:
     argList = '--startGPS %d --endGPS %d --IFO %s --fMin %d --fMax %d --SFTs %s/SFT*' % (analysisStartTime,analysisEndTime,ifo,thisStartFreq,thisEndFreq,pathToSFTs)
  tagStringOut = '%s_%i' % (tagString, nodeCount)  
  dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(specAvgJobName,argList,tagStringOut))
  if (createSFTs):  
    dagFID.write('PARENT %s CHILD %s\n'%(sftDAGSUBJobName,specAvgJobName))
  if (makeMatlabPlots):
     runMatlabPlotScriptJobName = 'runMatlabPlotScript_%i' % nodeCount
     dagFID.write('JOB %s runMatlabPlotScript.sub\n' % runMatlabPlotScriptJobName)
     dagFID.write('RETRY %s 3\n' % runMatlabPlotScriptJobName)
     inputFileName = 'spec_%d.00_%d.00_%s_%d_%d' % (thisStartFreq,thisEndFreq,ifo,analysisStartTime,analysisEndTime)
     # Matlab will append .pdf and .png to outputFileName to save plots:
     outputFileName = '%s/%s' % (plotOutputPath, inputFileName)     
     effTBase = timeBaseline/180.0
     effTBaseFull = timeBaseline
     deltaFTicks = 5
     taveFlag = 1
     if (htmlReferenceDir != None):
        referenceFileName = '%s/spec_%d.00_%d.00_%s' % (htmlReferenceDir,thisStartFreq,thisEndFreq,htmlRefIFOEpoch)
     else:
        referenceFileName = 'none'
     argList = '%s %s %s %d %d %d %d %d %d %s' % (inputFileName,outputFileName,channelName,effTBase,deltaFTicks,taveFlag,effTBaseFull,thresholdSNR,coincidenceDeltaF,referenceFileName)
     tagStringOut = '%s_%i' % (tagString, nodeCount)  
     dagFID.write('VARS %s argList="%s" tagstring="%s"\n'%(runMatlabPlotScriptJobName,argList,tagStringOut))
     dagFID.write('PARENT %s CHILD %s\n'%(specAvgJobName,runMatlabPlotScriptJobName))
  if (htmlFilename != None):
     inputFileName = 'spec_%d.00_%d.00_%s_%d_%d' % (thisStartFreq,thisEndFreq,ifo,analysisStartTime,analysisEndTime)
     # Matlab will append .pdf and .png to outputFileName to save plots.
     htmlFID.write('  <tr>\n')
     htmlFID.write('  <td style="vertical-align: top;">\n')
     htmlFID.write('    <a href="%s.pdf"><img alt="" src="%s.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n' % (inputFileName,inputFileName))
     htmlFID.write('  </td>\n')
     if (htmlReferenceDir != None):
        referenceFileName = '%s/spec_%d.00_%d.00_%s' % (htmlReferenceDir,thisStartFreq,thisEndFreq,htmlRefIFOEpoch)
        htmlFID.write('  <td style="vertical-align: top;">\n')
        htmlFID.write('    <a href="%s.pdf"><img alt="" src="%s.png" style="border: 0px solid ; width: 576px; height: 432px;"></a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('  </td>\n')     
     htmlFID.write('  </tr>\n')
     # Next row with links to text files:
     htmlFID.write('  <tr>\n')
     htmlFID.write('  <td style="vertical-align: top;">\n')
     htmlFID.write('    SFT Timestamps: <a href="%s_timestamps">%s_timestamps</a><br>\n' % (inputFileName,inputFileName))
     htmlFID.write('    Freq. vs Power: <a href="%s.txt">%s.txt</a><br>\n' % (inputFileName,inputFileName))
     htmlFID.write('    Freq. vs Power (Sorted): <a href="%s_sorted.txt">%s_sorted.txt</a><br>\n' % (inputFileName,inputFileName))
     if (htmlReferenceDir != None) and (thresholdSNR > 0):
        htmlFID.write('    Coincident Lines: <a href="%s_lines.txt">%s_lines.txt</a><br>\n' % (inputFileName,inputFileName))
        htmlLinesFID.write('<br>\n')
        htmlLinesFID.write('<object data="%s_lines.txt" type="text/plain" style="width: 620px; height: 450px"></object>' % inputFileName)
        htmlLinesFID.write('<br>\n')
     htmlFID.write('  </td>\n')
     if (htmlReferenceDir != None):
        htmlFID.write('  <td style="vertical-align: top;">\n')
        htmlFID.write('    SFT Timestamps: <a href="%s_timestamps">%s_timestamps</a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('    Freq. vs Power: <a href="%s.txt">%s.txt</a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('    Freq. vs Power (Sorted): <a href="%s_sorted.txt">%s_sorted.txt</a><br>\n' % (referenceFileName,referenceFileName))
        htmlFID.write('  </td>\n')
     htmlFID.write('  </tr>\n')
  thisStartFreq = thisStartFreq + freqSubBand
  thisEndFreq = thisStartFreq + freqSubBand
  nodeCount = nodeCount + 1

# Close the DAG file
dagFID.close

if (htmlFilename != None):
  htmlFID.write('</tbody>\n')
  htmlFID.write('</table>\n')
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
     htmlLinesFID.close

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
