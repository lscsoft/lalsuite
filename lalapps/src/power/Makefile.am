## Process this file with automake to produce Makefile.in
# $Id$
POWERSRCS = power.c Power.h 
if DATAFLOW
DATAFLOWPROGRAMS = lalapps_snglBurstHistogram lalapps_burca lalapps_bread
else
DATAFLOWPROGRAMS = 
endif
if FRAME
FRAMEPROGRAMS = lalapps_power
else
FRAMEPROGRAMS = 
endif
bin_PROGRAMS = $(DATAFLOWPROGRAMS) $(FRAMEPROGRAMS)
if FRAME
lalapps_power_SOURCES = $(POWERSRCS)
endif
if DATAFLOW
lalapps_snglBurstHistogram_SOURCES = snglBurstHistogram.c
lalapps_burca_SOURCES = burca.c
lalapps_bread_SOURCES = bread.c
endif
EXTRA_DIST = $(POWERSRCS) 
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la
