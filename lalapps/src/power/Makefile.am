## Process this file with automake to produce Makefile.in
# $Id$
POWERSRCS = power.c power.h 
if LALMETAIO
LALMETAIOPROGRAMS = lalapps_snglBurstHistogram lalapps_burca lalapps_bread lalapps_binj lalapps_binj_find
else
LALMETAIOPROGRAMS = 
endif
if LALFRAME
LALFRAMEPROGRAMS = lalapps_power
else
LALFRAMEPROGRAMS = 
endif
bin_PROGRAMS = $(LALMETAIOPROGRAMS) $(LALFRAMEPROGRAMS)
if LALFRAME
lalapps_power_SOURCES = $(POWERSRCS)
endif
if LALMETAIO
lalapps_snglBurstHistogram_SOURCES = snglBurstHistogram.c
lalapps_burca_SOURCES = burca.c
lalapps_bread_SOURCES = bread.c
lalapps_binj_SOURCES = binj.c
lalapps_binj_find_SOURCES = binj_find.c
endif
EXTRA_DIST = $(POWERSRCS) 
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
