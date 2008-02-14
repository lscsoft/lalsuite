## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
frameprogs = lalapps_MakeSFTs
else
frameprogs =
endif

bin_PROGRAMS = $(frameprogs) lalapps_ConvertToSFTv2 lalapps_ConvertChunkToSFTv2 lalapps_LISAmakeSFTs

if FRAME
lalapps_MakeSFTs_SOURCES = MakeSFTs.c
endif

lalapps_ConvertToSFTv2_SOURCES = ConvertToSFTv2.c

lalapps_ConvertChunkToSFTv2_SOURCES = ConvertChunkToSFTv2.c

lalapps_LISAmakeSFTs_SOURCES = LISAmakeSFTs.c readxml.c readxml.h ezxml.c ezxml.h
