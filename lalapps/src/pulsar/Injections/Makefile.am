## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

AM_LDFLAGS = -all-static

bin_PROGRAMS = makefakedata_v2 lalapps_Makefakedata compareSFTs dumpSFT makefakedatabinary makefakedata_test s3inject

lalapps_Makefakedata_SOURCES = makefakedata_v4.c	

makefakedata_v2_SOURCES =  makefakedata_v2.c

compareSFTs_SOURCES = compareSFTs.c

dumpSFT_SOURCES = dumpSFT.c

makefakedatabinary_SOURCES = makefakedatabinary.c

makefakedata_test_SOURCES = makefakedata_test.c

s3inject_SOURCES = s3inject.c

