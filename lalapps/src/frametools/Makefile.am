## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
FRAMEPROGRAMS = lalapps_animate lalapps_calfacs lalapps_computespec \
               lalapps_frinfo lalapps_frread lalapps_frextr lalapps_frjoin
#lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_animate_SOURCES = animate.c
lalapps_calfacs_SOURCES = calfacs.c
lalapps_computespec_SOURCES = computespec.c
lalapps_frinfo_SOURCES = info.c
lalapps_frread_SOURCES = read.c
lalapps_frextr_SOURCES = extr.c
lalapps_frjoin_SOURCES = join.c
else
FRAMEPROGRAMS =
endif
if LALFRAME
LALFRAMEPROGRAMS = lalapps_frview
lalapps_frview_SOURCES = view.c
else
LALFRAMEPROGRAMS =
endif
bin_PROGRAMS = $(FRAMEPROGRAMS) $(LALFRAMEPROGRAMS)
EXTRA_DIST = stripAdc.c lalDebugLevel.c animate.c calfacs.c computespec.c
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
