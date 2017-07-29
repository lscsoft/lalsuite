MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

AM_CPPFLAGS += -I$(top_srcdir)/src/pulsar/Fstatistic

bin_PROGRAMS = lalapps_SFTclean lalapps_HoughValidate \
	lalapps_DriveHoughMulti lalapps_MCInjectHoughMulti lalapps_MultiWeights \
	lalapps_ValidateHoughMulti lalapps_ValidateChi2Test lalapps_MCInjectHoughMultiChi2Test \
	lalapps_HoughValidateAM lalapps_HoughMismatch

lalapps_SFTclean_SOURCES = SFTclean.c
lalapps_HoughValidate_SOURCES = HoughValidate.c MCInjectComputeHough.h DriveHoughColor.h PeakSelect.h PeakSelect.c SFTfileIOv1.h SFTfileIOv1.c
lalapps_DriveHoughMulti_SOURCES = DriveHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h \
	 FstatToplist.c FstatToplist.h \
	 HeapToplist.c

lalapps_ValidateHoughMulti_SOURCES = ValidateHoughMulti.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
lalapps_MultiWeights_SOURCES = MultiWeights.c  DriveHoughColor.h
lalapps_MCInjectHoughMulti_SOURCES = MCInjectHoughMulti.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h
lalapps_ValidateChi2Test_SOURCES = ValidateHoughMultiChi2Test.c  DriveHoughColor.h PeakSelect.c PeakSelect.h MCInjectHoughMulti.h
lalapps_MCInjectHoughMultiChi2Test_SOURCES = MCInjectHoughMultiChi2Test.c MCInjectHoughMulti.h DriveHoughColor.h  PeakSelect.c PeakSelect.h

lalapps_HoughValidateAM_SOURCES = HoughValidateAM.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c SFTfileIOv1.h SFTfileIOv1.c
lalapps_HoughMismatch_SOURCES = HoughMismatch.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c SFTfileIOv1.h SFTfileIOv1.c
