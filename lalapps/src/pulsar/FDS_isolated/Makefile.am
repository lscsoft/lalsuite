
## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
AM_LDFLAGS = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
FRAMEPROGRAMS = lalapps_FrComputeFStatistic
else
FRAMEPROGRAMS = 
endif

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF lalapps_polka lalapps_FstatShapeTestLAL ComputeFStatisticBinary  compareFstats $(FRAMEPROGRAMS) ComputeFStatistic_v2 lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh

CFSsources = 	ComputeFStatistic.h		\
		rngmed.c rngmed.h clusters.c clusters.h  	\
		DopplerScan.c DopplerScan.h NeighboringGridPoints.c

lalapps_ComputeFStatistic_SOURCES = ComputeFStatistic.c $(CFSsources)
lalapps_ComputeFStatistic_LDFLAGS = -all-static

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

##----------------------------------------------------------------------
## some optional BOINC-related stuff follows here:
EXTRA_PROGRAMS = cfsBOINC cfsBOINCnographics

cfsBOINC_SOURCES = cfsBOINC.cc $(CFSsources) \
	screensaver/src/starlist.C screensaver/src/starsphere.C \
	screensaver/src/user_text.C screensaver/src/pulsar_list.C \
	screensaver/src/search_info.C
cfsBOINC_CPPFLAGS = $(AM_CPPFLAGS) -I$(HOME_OPENGL)/include -I$(HOME_BOINC)/include/BOINC -DUSE_BOINC 
cfsBOINC_LDFLAGS =  $(AM_LDFLAGS) -static \
	-L$(HOME_BOINC)/lib -lboinc -lboinc_api -lboinc_graphics_api \
	-L$(HOME_OPENGL)/lib -L/usr/X11R6/lib -lboincglut -lGL -lGLU -lXmu  

# non-graphic version
cfsBOINCnographics_SOURCES = cfsBOINC.cc $(CFSsources)
cfsBOINCnographics_CPPFLAGS = $(AM_CPPFLAGS) -I$(HOME_OPENGL)/include -I$(HOME_BOINC)/include/BOINC \
	-DUSE_BOINC -DNO_BOINC_GRAPHICS
cfsBOINCnographics_LDFLAGS = $(AM_LDFLAGS) -L$(HOME_BOINC)/lib -all-static -lboinc -lboinc_api 

## trick automake into using g++ for ComputeFStatistic.c (for boinc-includes to work):
cfsBOINC.cc: ComputeFStatistic.c
	ln -s ComputeFStatistic.c cfsBOINC.cc
