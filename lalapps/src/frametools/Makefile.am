## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
bin_PROGRAMS = lalapps_animate lalapps_calfacs lalapps_computespec
lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_animate_SOURCES = animate.c
lalapps_calfacs_SOURCES = calfacs.c
lalapps_computespec_SOURCES = computespec.c
else
bin_PROGRAMS = 
endif
EXTRA_DIST = stripAdc.c lalDebugLevel.c animate.c calfacs.c computespec.c
INCLUDES = -I$(top_srcdir)/src
