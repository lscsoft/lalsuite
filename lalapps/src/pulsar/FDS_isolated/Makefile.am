## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
AM_LDFLAGS = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
FRAMEPROGRAMS = lalapps_FrComputeFStatistic
else
FRAMEPROGRAMS = 
endif

bin_PROGRAMS = lalapps_ComputeFStatistic lalapps_FstatShapeTest lalapps_SemiAnalyticF \
	lalapps_polka lalapps_FstatShapeTestLAL compareFstats ComputeFStatistic_v2 \
	lalapps_makeInvetofile lalapps_extractSFTband lalapps_FindSh  \

EXTRA_PROGRAMS = cfsBOINC cfsBOINCnographics $(FRAMEPROGRAMS) 


CFSsources = 	ComputeFStatistic.h	\
		clusters.c clusters.h  	\
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

lalapps_FindSh_SOURCES=FindSh.c FindSh.h clusters.c clusters.h

compareFstats_SOURCES = compareFstats.c

ComputeFStatistic_v2_SOURCES = ComputeFStatistic_v2.c clusters.c clusters.h  \
	DopplerScan.c DopplerScan.h NeighboringGridPoints.c

##----------------------------------------------------------------------
## some optional BOINC-related stuff follows here:
HOME_X11 = /usr/X11R6/lib
cfsBOINC_SOURCES = cfsBOINC.cc $(CFSsources) \
	screensaver/src/starlist.C screensaver/src/starsphere.C \
	screensaver/src/user_text.C screensaver/src/pulsar_list.C \
	screensaver/src/search_info.C
cfsBOINC_CPPFLAGS = -Wall $(AM_CPPFLAGS) -I$(HOME_OPENGL)/include -I$(HOME_BOINC)/include/BOINC -DUSE_BOINC -D_REENTRANT
cfsBOINC_LDFLAGS =  $(AM_LDFLAGS) -static \
	-L$(HOME_BOINC)/lib -lboinc -lboinc_api -lboinc_graphics_api \
	-L$(HOME_OPENGL)/lib -L$(HOME_X11) -lboincglut -lpthread -lGL -lGLU -lXmu -lX11 -lXext -ldl -lm -lpthread -lXi -lXmu -lX11

#cfsBOINC_LDFLAGS = -static -L$(HOME_BOINC)/lib -L$(HOME_OPENGL)/lib -lboinc \
#-print-file-name=libpthread.so -lboinc_api -lboinc_graphics_api -lGL -lboincglut -L/usr/X11R6/lib -lGL -lGLU -lXmu

cfsBOINCnographics_SOURCES = cfsBOINC.cc $(CFSsources)
cfsBOINCnographics_CPPFLAGS = $(AM_CPPFLAGS) -I$(HOME_OPENGL)/include -I$(HOME_BOINC)/include/BOINC -DUSE_BOINC -DNO_BOINC_GRAPHICS -D_REENTRANT
cfsBOINCnographics_LDFLAGS = $(AM_LDFLAGS) -L$(HOME_BOINC)/lib -static -lboinc -lboinc_api 

## trick make into using c++ for ComputeFStatistic.c (for boinc-includes):
cfsBOINC.cc: ComputeFStatistic.c
	ln -s ComputeFStatistic.c cfsBOINC.cc
