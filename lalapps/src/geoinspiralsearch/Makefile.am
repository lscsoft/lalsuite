## Process this file with automake to produce Makefile.in
# $Id$
bin_PROGRAMS = lalapps_GEO600InspiralSearch lalapps_GEO600InspiralSearchSnglProc
#bin_PROGRAMS = lalapps_GEO600InspiralSearch lalapps_GEO600InspiralSearchSnglProc lalapps_BankEfficiency
lalapps_GEO600InspiralSearch_SOURCES = GEO600InspiralSearch.c GEO600InspiralSearchFuncs.c GEO600InspiralSearchMasterDB.c GEO600InspiralSearchSlave.c GEO600FrameChannels.c 
lalapps_GEO600InspiralSearchSnglProc_SOURCES = GEO600InspiralSearchSnglProc.c GEO600InspiralSearchFuncs.c GEO600FrameChannels.c 
#lalapps_BankEfficiency_SOURCES = InspiralBankEfficiency.c InspiralBankEfficiencyMaster.c InspiralBankEfficiencySlave.c GEO600InspiralSearchFuncs.c GEO600FrameChannels.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

EXTRA_DIST = GEO600FrameChannels.h
