#!/usr/bin/env tclshexe
package require LDASJob
package require segments

# List of arguments:
# argv[1]: LDAS username (optional, use whoami if not provided)
# argv[2]: block number to start (optional) 

# get user name
if { $argc < 1 } {
    set userName $env(USER)
} else {
    set userName [ lindex $argv 0 ]
}

#######################################################
#######################################################
# User defined parameters:
#######################################################
#######################################################
# LDAS system to use
set manager cit

# Comma separated list of pre-processing filters.
# If more than one filter is given, they are applied serially.
# For each filter in the list, two files named filter_a.ilwd and 
# filter_b.ilwd must be present in the current directory. Otherwise,
# will try to get these files from $filtroot on the LDAS system.
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

# type of segments: must be a file which can be read
# by ligotools segment library
set dataType "S2H1H2L1v03_segs.txt"

# channel to process (can be a comma separated list, in which
# case the code is run on each channel in a loop)
# use channel in the ETG parameters
set channel "H1:LSC-AS_Q"

#######################################################
# Injection data
#######################################################

# Comma separated list of injection waveforms (can be empty)
# For each waveform in the list, at least on file of the form waveform_p.ilwd
# or waveform_c.ilwd must be present in the current directory.
# If one of the two files is missing, the file Zeros.ilwd is used in its
# place, and must also be present.
# Use the matlab file lalapps/src/burstwrapper/generateSGilwd.m for an 
# example of how to generate these files.
# All waveforms are given in strain units.
set waveforms ""

# NOTE on parameters: they can be
# o single values
# o list enclosed in parantheses: (1,2,3,4)
# o linear range: (\\\\\[min,max,number\\\\\]), ex: (\\\\\[1,6,5\\\\\])
# o log range: (\\\\\[min,max,number,l\\\\\]), ex: (\\\\\[1e-19,1e-16,9,l\\\\\])
# o random range: (\\\\\[min,max,number,r\\\\\]), ex: (\\\\\[32768,5865472,10,r\\\\\])
# The many escape characters are necessary to deal with the
# many levels of TCL code. Isn't LDAS lovely?
# o a vector in an ilwd file: _filename_, where filename.ilwd is passed to the wrapper. 
# o a matrix in an ilwd file: __filename__, where filename.ilwd is passed to the wrapper. See lalapps/src/burstwrapper/ShellInjection.c for an example of how this is generated.

# injection strain amplitude.
# Numerical factor that multiplies the waveform when it is added to the
# ADC data after its filtering by respfilt in the datacondAPI, in order
# to convert it from strain to ADC.
set injAmp ""

# injection alpha (if outside [0,2pi], inject at zenith)
# Follows LAL convention.
set injAlpha ""

# injection delta (if outside [0,2pi], inject at zenith)
# Follows LAL convention.
set injDelta ""

# injection polarization angle psi
# Follows LAL convention.
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

# Parameters of the ETG
# In this form, the parameters are passed from the file H1etg.ilwd.
# This file was produced using lalapps/src/burstwrapper/ilwd_burstdso.m 
set ETGParameters "__H1etg__"

# As an alternative, each ETG parameter can be explicitely given:
#set TFCchn "channel"
#set TFCThr "0.023357"
#set TFCWin "1"
#set TFCThrMethod "0"
#set TFCdata "0"
#set TFCTRez "0.125"
#set TFCFmin "-130.0"
#set TFCFmax "-400.0"
#set TFCalpha "-1.0"
#set TFCsigma "2"
#set TFCdelta "0,0,0,0,0,0,0,0,0,0,0,0,0,0,0"

#######################################################
# Run parameters
#######################################################

# number of nodes
set NNodes 4

# frame type
set FType RDS_R_L1

# Burst MDC frames to add to the raw data in the datacondAPI
set MDCFrames SG10

# sampling frequency (integer) [Hz]
set samplingF 16384

# maximum number of times to retry a failed job
set maxRetry 2

# set to 1 to skip the burst parameter estimation function
set BurstOutput 0

# do 1 segment in ...
set do1inN 1

# output type (0:XML,1:binXML,2:bin to file)
# 0 and 1 not well tested.
set OutputType 2

# if OutputType == 2, path to file
# if "LOCAL:/usr1/...", save to local disk on node running the job. 
# Be careful with sending data to a shared filesystem, as severe
# reductions in performance might occur.
# The following is an example of using a shared filesystem:
set PrebinFile "/dso-test/LDAS_Driver_Example_$userName/"

# If you use "LOCAL:", you must submit from a machine that sees all the
# nodes. The script will try to create directories /data/node*/<dir>,
# for PrebinFile set to "LOCAL:/usr1/<dir>".
#set PrebinFile "LOCAL:/usr1/$userName/LDAS_Driver_Example_$userName/"

# if OutputType == 2, set to 1 to gzip the files
set Zip 0

# if OutputType == 2, set to 1 to move the files to current directory
set moveIt 0

#######################################################
#######################################################
# End of user defined parameters.
#######################################################
#######################################################


# This calls the main routines that submit the jobs to the
# requested LDAS system.
source $env(BURSTWRAPPER)/BurstDriverEngineV3MDC.tcl
