## Process this file with automake to produce Makefile.in
# $Id$
ZMDFTSRCS = zmdft.c zmdft.h 
bin_PROGRAMS = lalapps_zmdft lalapps_build_zm_frames
lalapps_zmdft_SOURCES = $(ZMDFTSRCS)
lalapps_build_zm_frames_SOURCES = build_zm_frames.c
EXTRA_DIST = $(ZMDFTSRCS) 
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
