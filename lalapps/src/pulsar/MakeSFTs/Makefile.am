## Process this file with automake to produce Makefile.in
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

if FRAME
frameprogs = lalapps_MakeSFTs lalapps_MakeCleanSFTs lalapps_MakeCleanSFTsNObil2 lalapps_MakePSSCleanSFTs
else
frameprogs =
endif

bin_PROGRAMS = $(frameprogs) lalapps_ConvertToSFTv2 lalapps_LISAmakeSFTs

if FRAME
lalapps_MakeSFTs_SOURCES = MakeSFTs.c
lalapps_MakeCleanSFTs_SOURCES = MakeCleanSFTs.c
lalapps_MakeCleanSFTsNObil2_SOURCES=MakeCleanSFTsNObil2.c
lalapps_MakePSSCleanSFTs_SOURCES = MakePSSCleanSFTs.c XLALPSSInterface.c
lalapps_MakePSSCleanSFTs_LDADD = -lsfdb -lpss -lastro -lnovas
endif

lalapps_ConvertToSFTv2_SOURCES = ConvertToSFTv2.c

lalapps_LISAmakeSFTs_SOURCES = LISAmakeSFTs.c readxml.c readxml.h ezxml.c ezxml.h
