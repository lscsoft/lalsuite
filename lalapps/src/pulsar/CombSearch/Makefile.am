## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -I$(top_srcdir)/misc -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la


## make sure LALAPPSgitID.h is re-generated at every 'make'
#.PHONY: gitID
#gitID:
#	cd $(top_builddir)/misc && $(MAKE) gitID
#BUILT_SOURCES = gitID

bin_PROGRAMS = lalapps_CombSearch

lalapps_CombSearch_SOURCES = CombSearch.c


## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common

