## Process this file with automake to produce Makefile.in
# $Id$

# for using lalapps-defines + functions
INCLUDES = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/src/pulsar/FDS_isolated -I$(top_srcdir)/misc
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

## make sure LALAPPSgitID.h is re-generated at every 'make'
.PHONY: gitID
gitID:
	cd $(top_builddir)/misc && $(MAKE) gitID
BUILT_SOURCES = gitID

bin_PROGRAMS = lalapps_HierarchSearchGCT

lalapps_HierarchSearchGCT_SOURCES = HierarchSearchGCT.c HierarchSearchGCT.h \
                          GCTtoplist.c GCTtoplist.h \
                          ../FDS_isolated/HeapToplist.c ../FDS_isolated/HeapToplist.h \
                          ../FDS_isolated/ComputeFStatistic.h ../FDS_isolated/Fstat_v3.c ComputeFstat_RS.c
