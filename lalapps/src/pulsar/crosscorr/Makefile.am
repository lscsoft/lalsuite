## Process this file with automake to produce Makefile.in

# for using lalapps-defines + functions
AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps  -I$(top_srcdir)/src/pulsar/FDS_isolated -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

bin_PROGRAMS = lalapps_pulsar_crosscorr lalapps_pulsar_crosscorr_v2

lalapps_pulsar_crosscorr_SOURCES = pulsar_crosscorr.c pulsar_crosscorr.h

lalapps_pulsar_crosscorr_v2_SOURCES = pulsar_crosscorr_v2.c 


## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/gnuscripts/vcsID.common
