## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
bin_PROGRAMS = lalapps_stripAdc lalapps_animate
lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
else
bin_PROGRAMS = lalapps_animate
endif
lalapps_animate_SOURCES = animate.c
EXTRA_DIST = stripAdc.c lalDebugLevel.c
INCLUDES = -I$(top_srcdir)/src
