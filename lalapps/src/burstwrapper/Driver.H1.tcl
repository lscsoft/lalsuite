#!/usr/bin/env tclshexe
package require LDASJob
package require segments

# List of arguments:
# argv[1]: LDAS username
# argv[2]: block number to start (optional) 

set userName [ lindex $argv 0 ]

#######################################################
#######################################################
# User defined parameters:
#######################################################
#######################################################
# LDAS system to use
set manager cit

# Comma separated list of pre-processing filters.
# For each filter in the list, two files named filter_a.ilwd and 
# filter_b.ilwd must be present in the current directory. Otherwise,
# will try to get these files from filtroot on the LDAS system.
# Use \$bId for the block number, \$PbId for Peter Shawhan's segment number
# Use \$IFO2 for IFO name
set prefilters "S2_HPF_6_100"

# FIR whitening filters (with b coefficients only)
set prefiltersB "S2_LPEF_\${IFO2}_\${PbId}_4_1_1024"

# path of ldas filter cache (if desired)
set filtroot "/ldas_outgoing/jobs/ilwd/filters/S2"

# duration of analyzed segment (s)
set duration 300

# amount of time to skip at beginning of segment after filtering (s)
# (duration+1) seconds will be acquired from beginning of segment, will
# be filtered, and the TransientSkip first seconds will be skipped.
set TransientSkip 1

# type of segments: 0=playground, 1=production, else=use that file
set dataType "S2H1H2L1v03_segs.txt"

# channel to process (can be a comma separated list)
# use channel in the ETG parameters
set channel "H1:LSC-AS_Q"

#######################################################
# Injection data
#######################################################.txt"

# channel to process (can be a comma separated list)
# use channel in the ETG parameters
set channel "H1:LSC-AS_Q"

#######################################################
# Injection data
#######################################################

# Comma separated list of injection waveforms (can be empty)
# For each waveform in the list, at least on file of the form waveform_p.ilwd
# or waveform_c.ilwd must be present in the current directory.
# If one of the two files is missing, the file Zeros.ilwd is used in its
# place.
set waveforms ""

# NOTE on parameters: they can be
# o single values
# o list enclosed in parantheses
# o linear range: (\\\\\[1,6,5\\\\\])
# o log range: (\\\\\[1e-19,1e-16,9,l\\\\\])
# o random range: (\\\\\[32768,5865472,10,r\\\\\])
# The many escape characters are necessary to deal with the
# many levels of TCL code. Isn't LDAS lovely?
# o a vector in an ilwd file: _filename_, where filename.ilwd is passed to the wrapper
# o a matrix in an ilwd file: __filename__, where filename.ilwd is passed to the wrapper

# injection amplitude
set injAmp ""

# injection alpha (if outside [0,2pi], inject at zenith)
set injAlpha ""

# injection delta (if outside [0,2pi], inject at zenith)
set injDelta ""

# injection polarization angle psi
set injPsi ""

# total number of injections to perform per segment
set injN 0

# times for injections to be done at one time in one segment (samples)
set injTimes ""

#######################################################
# ETG parameters
#######################################################

# name of the ETG to use
set ETG "TFCLUSTERS"

#set ETGParameters "__H1etg__"

set TFCchn "channel"
set TFCThr "0.023357"
set TFCWin "1"
set TFCThrMethod "0"
set TFCdata "0"
set TFCTRez "0.125"
set TFCFmin "-130.0"
set TFCFmax "-400.0"
set TFCalpha "-1.0"
set TFCsigma "2"
set TFCdelta "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"

#######################################################
# Run parameters
#######################################################

# number of nodes
set NNodes 4

# frame type
set FType RDS_R_L1

# MDC frames to add
set MDCFrames SG10

# sampling frequency (integer)
set samplingF 16384

# maximum number of times to retry a failed job
set maxRetry 2

# set to 1 to skip the burst parameter estimation function
set BurstOutput 0

# do 1 segment in ...
set do1inN 1

# output type (0:XML,1:binXML,2:bin to file)
set OutputType 2

# if OutputType == 2, path to file
# if "LOCAL:/usr1/...", save to local disk
set PrebinFile "LOCAL:/usr1/jsylvest/TEST_$userName/"

# if OutputType == 2, set to 1 to gzip the files
set Zip 0

# if OutputType == 2, set to 1 to move the files to current directory
set moveIt 0

#######################################################
#######################################################
# End of user defined parameters.
#######################################################
#######################################################


source $env(BURSTWRAPPER)/BurstDriverEngineSA.tcl
