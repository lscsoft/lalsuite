## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/src/pulsar/FDS_isolated
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_HoughValidateAM lalapps_HoughMismatch lalapps_HierarchicalSearch

EXTRA_PROGRAMS = lalapps_HierarchicalSearch_v2 lalapps_HierarchicalSearchGC

lalapps_HoughValidateAM_SOURCES = HoughValidateAM.c ../src/MCInjectHoughS2.h \
	                  ../src/SFTbin.h ../src/SFTbin.c \
			  ../src/DriveHoughColor.h ../src/PeakSelect.h \
			  ../src/PeakSelect.c

lalapps_HoughMismatch_SOURCES = HoughMismatch.c ../src/MCInjectHoughS2.h \
			../src/SFTbin.h ../src/SFTbin.c \
			../src/DriveHoughColor.h ../src/PeakSelect.h \
			../src/PeakSelect.c

lalapps_HierarchicalSearch_SOURCES = HierarchicalSearch.c HierarchicalSearch.h \
			  StackSlideFstat.c StackSlideFstat.h \
			  HoughFStatToplist.c HoughFStatToplist.h \
			  ../../FDS_isolated/HeapToplist.c ../../FDS_isolated/HeapToplist.h \
			  ../../FDS_isolated/ComputeFStatistic.h

lalapps_HierarchicalSearch_v2_SOURCES = HierarchicalSearch_v2.c HierarchicalSearch.h StackSlideFstat.c StackSlideFstat.h \
			  ../../FDS_isolated/FstatToplist.c ../../FDS_isolated/FstatToplist.h \
			  ../../FDS_isolated/HeapToplist.c ../../FDS_isolated/HeapToplist.h \
			  ../../FDS_isolated/ComputeFStatistic.h

lalapps_HierarchicalSearchGC_SOURCES = HierarchicalSearchGC.c HierarchicalSearchGC.h \
			  StackSlideFstat.c StackSlideFstat.h \
			  HoughFStatToplist.c HoughFStatToplist.h \
			  ../../FDS_isolated/HeapToplist.c ../../FDS_isolated/HeapToplist.h \
			  ../../FDS_isolated/ComputeFStatistic.h
