## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
#INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
#LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## not usable before automake-1.6:
## AM_LDFLAGS = -all-static

bin_PROGRAMS = HoughValidateAM

HoughValidateAM_SOURCES = HoughValidateAM.c HoughValidateAM.h ../src/SFTbin.h ../src/SFTbin.c ../src/DriveHoughColor.h ../src/PeakSelect.h ../src/PeakSelect.c
