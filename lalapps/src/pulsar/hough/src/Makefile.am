## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## not usable before automake-1.6:
## AM_LDFLAGS = -all-static

bin_PROGRAMS = NstarClean SFTclean MCInjectS2 DriveHough_v3 HoughValidate ComputePSD DriveHoughMulti MCInjectHoughMulti MultiWeights ValidateHoughMulti

NstarClean_SOURCES = NstarClean.c
SFTclean_SOURCES = SFTclean.c 
ComputePSD_SOURCES = ComputePSD.c 
MCInjectS2_SOURCES = MCInjectHoughS2.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
DriveHough_v3_SOURCES = DriveHough_v3.c  DriveHoughColor.h PeakSelect.c PeakSelect.h
HoughValidate_SOURCES = HoughValidate.c MCInjectComputeHough.h DriveHoughColor.h PeakSelect.h PeakSelect.c
DriveHoughMulti_SOURCES = DriveHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h
ValidateHoughMulti_SOURCES = ValidateHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
MultiWeights_SOURCES = MultiWeights.c  DriveHoughColor.h
MCInjectHoughMulti_SOURCES = MCInjectHoughMulti.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h
