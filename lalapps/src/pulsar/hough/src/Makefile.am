## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
#INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
#LDADD = $(top_builddir)/src/lalapps/liblalapps.la

AM_LDFLAGS = -all-static

bin_PROGRAMS = NstarClean SFTclean MCInjectS2 DriveHough_v3 HoughValidate 

NstarClean_SOURCES = NstarClean.c SFTbin.h SFTbin.c
SFTclean_SOURCES = SFTclean.c SFTbin.h SFTbin.c
MCInjectS2_SOURCES = MCInjectHoughS2.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
DriveHough_v3_SOURCES = DriveHough_v3.c  DriveHoughColor.h SFTbin.c PeakSelect.c
HoughValidate_SOURCES = HoughValidate.c MCInjectComputeHough.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
