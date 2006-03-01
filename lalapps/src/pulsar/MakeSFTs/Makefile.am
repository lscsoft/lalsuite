## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
	frameprogs = lalapps_MakeSFTs
	lalapps_MakeSFTs_SOURCES = MakeSFTs.c
else
	frameprogs = 
endif

bin_PROGRAMS = $(frameprogs) ConvertToSFTv2

ConvertToSFTv2_SOURCES = ConvertToSFTv2.c
