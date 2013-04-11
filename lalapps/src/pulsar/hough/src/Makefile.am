## Process this file with automake to produce Makefile.in

# for using lalapps-defines + functions
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(top_srcdir)/src/pulsar/FDS_isolated
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_NstarClean lalapps_SFTclean lalapps_MCInjectS2 lalapps_DriveHough_v3 lalapps_HoughValidate \
	lalapps_DriveHoughMulti lalapps_MCInjectHoughMulti lalapps_MultiWeights \
	lalapps_ValidateHoughMulti lalapps_ValidateChi2Test lalapps_MCInjectHoughMultiChi2Test

lalapps_NstarClean_SOURCES = NstarClean.c
lalapps_SFTclean_SOURCES = SFTclean.c
lalapps_MCInjectS2_SOURCES = MCInjectHoughS2.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
lalapps_DriveHough_v3_SOURCES = DriveHough_v3.c  DriveHoughColor.h PeakSelect.c PeakSelect.h
lalapps_HoughValidate_SOURCES = HoughValidate.c MCInjectComputeHough.h DriveHoughColor.h PeakSelect.h PeakSelect.c
lalapps_DriveHoughMulti_SOURCES = DriveHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h \
	 ../../FDS_isolated/FstatToplist.c ../../FDS_isolated/FstatToplist.h \
	 ../../FDS_isolated/HeapToplist.c ../../FDS_isolated/HeapToplist.h

lalapps_ValidateHoughMulti_SOURCES = ValidateHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
lalapps_MultiWeights_SOURCES = MultiWeights.c  DriveHoughColor.h
lalapps_MCInjectHoughMulti_SOURCES = MCInjectHoughMulti.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h
lalapps_ValidateChi2Test_SOURCES = ValidateHoughMultiChi2Test.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
lalapps_MCInjectHoughMultiChi2Test_SOURCES = MCInjectHoughMultiChi2Test.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h


## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
