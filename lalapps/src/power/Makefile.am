## Process this file with automake to produce Makefile.in
# $Id$
POWERSRCS = power.c \
  Power.h 
if FRAME
bin_PROGRAMS = lalapps_power
lalapps_power_SOURCES = $(POWERSRCS)
else
bin_PROGRAMS =
endif
EXTRA_DIST = $(POWERSRCS) 
INCLUDES = -I$(top_srcdir)/src
