## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps -DPKG_DATA_DIR='"$(pkgdatadir)/"' @AM_CPPFLAGS@
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_TEMPOcomparison

lalapps_TEMPOcomparison_SOURCES = TEMPOcomparison.c

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
