## Process this file with automake to produce Makefile.in
# $Id$
ZMDFTSRCS = zmdft.c zmdft.h 
ZMCOMPARESRCS = zmcompare.c
if LALFRAME
LALFRAMEPROGRAMS = lalapps_zmdft lalapps_build_zm_frames
else
LALFRAMEPROGRAMS =
endif
bin_PROGRAMS = $(LALFRAMEPROGRAMS) lalapps_zmcompare
if LALFRAME
lalapps_zmdft_SOURCES = $(ZMDFTSRCS)
lalapps_zmcompare_SOURCES = $(ZMCOMPARESRCS)
lalapps_build_zm_frames_SOURCES = build_zm_frames.c
endif
EXTRA_DIST = $(ZMDFTSRCS) $(ZMCOMPARESRCS)
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
