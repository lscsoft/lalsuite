## Process this file with automake to produce Makefile.in
# $Id$
POWERSRCS = power.c InitSearch.c FinalizeSearch.c ConditionData.c \
  Power.h InitSearch.h FinalizeSearch.h ConditionData.h ApplySearch.h \
  LALWrapperConfig.h LALWrapperInterface.h \
  wrapperInterfaceDatatypes.h wrapperInterface.h
if FRAME
bin_PROGRAMS = lalapps_power
lalapps_power_SOURCES = $(POWERSRCS)
else
bin_PROGRAMS =
endif
EXTRA_DIST = $(POWERSRCS)
INCLUDES = -I$(top_srcdir)/src
