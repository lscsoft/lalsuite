## Process this file with automake to produce Makefile.in
# $Id$
if FRAME
bin_PROGRAMS = lalapps_power
lalapps_power_SOURCES = power.c InitSearch.c FinalizeSearch.c ConditionData.c
else
bin_PROGRAMS =
endif
EXTRA_DIST = power.c InitSearch.c FinalizeSearch.c ConditionData.c
INCLUDES = -I$(top_srcdir)/src
