## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
bin_PROGRAMS = lalapps_stripAdc lalapps_animate
lalapps_stripAdc_SOURCES = stripAdc.c lalDebugLevel.c
lalapps_animate_SOURCES = animate.c
else
bin_PROGRAMS = 
endif
EXTRA_DIST = stripAdc.c lalDebugLevel.c animate.c
INCLUDES = -I$(top_srcdir)/src
