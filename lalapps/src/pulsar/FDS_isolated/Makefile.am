## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## build everything completely static
##AM_LDFLAGS = -all-static

if FRAME
#FRAMEPROGRAMS = lalapps_FrComputeFStatistic
FRAMEPROGRAMS = 
else
FRAMEPROGRAMS = 
endif

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF lalapps_polka lalapps_FstatShapeTestLAL ComputeFStatisticBinary  compareFstats $(FRAMEPROGRAMS) ComputeFStatistic_v2 lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh

CFSsources = 	ComputeFStatistic.c ComputeFStatistic.h		\
		rngmed.c rngmed.h clusters.c clusters.h  	\
		DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_ComputeFStatistic_SOURCES = $(CFSsources)

lalapps_FrComputeFStatistic_SOURCES = FrComputeFStatistic.c FrComputeFStatistic.h	\
	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_FstatShapeTest_SOURCES = FstatShapeTest.c 
lalapps_FstatShapeTestLAL_SOURCES = FstatShapeTestLAL.c

lalapps_SemiAnalyticF_SOURCES =	SemiAnalyticF.c SemiAnalyticF.h

lalapps_polka_SOURCES = polka.c

lalapps_makeInvetofile_SOURCES=makeInvetofile.c

lalapps_extractSFTband_SOURCES=extractSFTband.c

lalapps_FindSh_SOURCES=FindSh.c FindSh.h clusters.c clusters.h rngmed.c rngmed.h

ComputeFStatisticBinary_SOURCES = ComputeFStatisticBinary.c ComputeFStatisticBinary.h rngmed.c rngmed.h clusters.c clusters.h

compareFstats_SOURCES = compareFstats.c

ComputeFStatistic_v2_SOURCES = ComputeFStatistic_v2.c 	rngmed.c rngmed.h clusters.c clusters.h  			\
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

## some optional BOINC-related stuff follows here:
EXTRA_PROGRAMS = cfsBOINC
cfsBOINC_SOURCES = $(CFSsources) screensaver/src/starlist.C screensaver/src/starsphere.C screensaver/src/user_text.C
cfsBOINC_LDFLAGS =  -L$(HOME_BOINC)/lib -L$(HOME_OPENGL)/lib -static -lboinc -lboinc_api -lboinc_graphics_api -lGL -lboincglut -lglut
cfsBOINC_CPPFLAGS = -I$(HOME_OPENGL)/include -I$(HOME_BOINC)/include -DUSE_BOINC ##-DNO_BOINC_GRAPHICS
