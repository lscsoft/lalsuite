## Process this file with automake to produce Makefile.in

LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_TEMPOcomparison

lalapps_TEMPOcomparison_SOURCES = TEMPOcomparison.c

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
