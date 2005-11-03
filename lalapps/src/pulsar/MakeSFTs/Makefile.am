## Process this file with automake to produce Makefile.in
# $Id$

INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if FRAME
  bin_PROGRAMS = lalapps_MakeSFTs
  lalapps_make_sfts_SOURCES = MakeSFTs.c
endif
