## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## not usable before automake-1.6:
## AM_LDFLAGS = -all-static

bin_PROGRAMS = HoughValidateAM HoughMismatch HierarchicalSearch

HoughValidateAM_SOURCES = HoughValidateAM.c ../src/MCInjectHoughS2.h \
	                  ../src/SFTbin.h ../src/SFTbin.c \
			  ../src/DriveHoughColor.h ../src/PeakSelect.h \
				  ../src/PeakSelect.c
HoughMismatch_SOURCES = HoughMismatch.c ../src/MCInjectHoughS2.h \
			../src/SFTbin.h ../src/SFTbin.c \
			../src/DriveHoughColor.h ../src/PeakSelect.h \
			../src/PeakSelect.c

HierarchicalSearch_SOURCES = HierarchicalSearch.c HierarchicalSearch.h \
			  ../../FDS_isolated/DopplerScan.c ../../FDS_isolated/DopplerScan.h \
			  ../../FDS_isolated/LogPrintf.c ../../FDS_isolated/LogPrintf.h	\
			  ../../FDS_isolated/FstatToplist.c ../../FDS_isolated/FstatToplist.h \
			  ../../FDS_isolated/ComputeFStatistic.h

