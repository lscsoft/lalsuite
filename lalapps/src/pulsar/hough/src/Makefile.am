## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
#INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
#LDADD = $(top_builddir)/src/lalapps/liblalapps.la

AM_LDFLAGS = -all-static

bin_PROGRAMS = SFTclean MCInject MCInjectS2 DriveHoughColor_v2 MCValidate 

SFTclean_SOURCES = SFTclean.c SFTbin.h SFTbin.c
MCInject_SOURCES = MCInjectComputeHoughII.c MCInjectComputeHough.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
MCInjectS2_SOURCES = MCInjectHoughS2.c MCInjectHoughS2.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
DriveHoughColor_v2_SOURCES = DriveHoughColor_v2.c  DriveHoughColor.h SFTbin.c PeakSelect.c
MCValidate_SOURCES = MCInjectValidate.c MCInjectComputeHough.h SFTbin.h SFTbin.c DriveHoughColor.h PeakSelect.h PeakSelect.c
