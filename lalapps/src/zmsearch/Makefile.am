## Process this file with automake to produce Makefile.in
# $Id$
ZMDFTSRCS = zmdft.c zmdft.h
ZMCOMPARESRCS = zmcompare.c
BUILDZMSRCS = build_zm_frames.c

if LALFRAME
if LALMETAIO
if LALBURST
LALPROGRAMS = lalapps_zmdft lalapps_build_zm_frames \
                   lalapps_zmcompare
lalapps_zmdft_SOURCES = $(ZMDFTSRCS)
lalapps_zmcompare_SOURCES = $(ZMCOMPARESRCS)
lalapps_build_zm_frames_SOURCES = ${BUILDZMSRCS}
endif
endif
endif
bin_PROGRAMS = $(LALPROGRAMS)

EXTRA_DIST = $(ZMDFTSRCS) $(ZMCOMPARESRCS) $(BUILDZMSRCS)
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
