## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/src/pulsar/FDS_isolated
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## not usable before automake-1.6:
## AM_LDFLAGS = -all-static

bin_PROGRAMS = NstarClean SFTclean MCInjectS2 DriveHough_v3 HoughValidate DriveHoughMulti MCInjectHoughMulti MultiWeights ValidateHoughMulti ValidateChi2Test MCInjectHoughMultiChi2Test

NstarClean_SOURCES = NstarClean.c
SFTclean_SOURCES = SFTclean.c 
MCInjectS2_SOURCES = MCInjectHoughS2.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
DriveHough_v3_SOURCES = DriveHough_v3.c  DriveHoughColor.h PeakSelect.c PeakSelect.h
HoughValidate_SOURCES = HoughValidate.c MCInjectComputeHough.h DriveHoughColor.h PeakSelect.h PeakSelect.c
DriveHoughMulti_SOURCES = DriveHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h \
	 ../../FDS_isolated/FstatToplist.c ../../FDS_isolated/FstatToplist.h \
	 ../../FDS_isolated/HeapToplist.c ../../FDS_isolated/HeapToplist.h

ValidateHoughMulti_SOURCES = ValidateHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
MultiWeights_SOURCES = MultiWeights.c  DriveHoughColor.h
MCInjectHoughMulti_SOURCES = MCInjectHoughMulti.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h
ValidateChi2Test_SOURCES = ValidateHoughMultiChi2Test.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
MCInjectHoughMultiChi2Test_SOURCES = MCInjectHoughMultiChi2Test.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h
