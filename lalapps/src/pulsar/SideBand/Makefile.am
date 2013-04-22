## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_GenerateSideBandTemplate lalapps_SideBandMCMC \
		lalapps_GetTimeStamps

lalapps_GenerateSideBandTemplate_SOURCES = GenerateSideBandTemplate.c SideBandUtils.c SideBand.h

lalapps_SideBandMCMC_SOURCES = SideBandMCMC.c SideBandUtils.c SideBand.h

lalapps_GetTimeStamps_SOURCES = GetTimeStamps.c

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
