## Process this file with automake to produce Makefile.in
# $Id$
ZMDFTSRCS = zmdft.c zmdft.h 
ZMCOMPARESRCS = zmcompare.c
BUILDZMSRCS = build_zm_frames.c

if LALFRAME
if LALMETAIO
LALFRAMEPROGRAMS = lalapps_zmdft lalapps_build_zm_frames \
                   lalapps_zmcompare
lalapps_zmdft_SOURCES = $(ZMDFTSRCS)
lalapps_zmcompare_SOURCES = $(ZMCOMPARESRCS)
lalapps_build_zm_frames_SOURCES = ${BUILDZMSRCS}
else
LALFRAMEPROGRAMS =
endif
else
LALFRAMEPROGRAMS =
endif
bin_PROGRAMS = $(LALFRAMEPROGRAMS)

EXTRA_DIST = $(ZMDFTSRCS) $(ZMCOMPARESRCS)
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
