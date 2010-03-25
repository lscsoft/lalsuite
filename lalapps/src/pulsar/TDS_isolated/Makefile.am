# Makefile for TDS_isolated
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALFRAME
  bin_PROGRAMS = lalapps_pulsar_parameter_estimation lalapps_heterodyne_pulsar
  TESTS = test_heterodyne_pulsar.sh
else
  bin_PROGRAMS =
  TESTS =
endif

lalapps_pulsar_parameter_estimation_SOURCES = pulsar_parameter_estimation.c pulsar_parameter_estimation.h

lalapps_heterodyne_pulsar_SOURCES = heterodyne_pulsar.c heterodyne_pulsar.h

EXTRA_DIST = test_heterodyne_pulsar.sh make_frame_cache H1response.txt H-CW_Injection-875206560-120.gwf

TESTS_ENVIRONMENT = \
LALSUITE_TOP_SRCDIR=$(LALSUITE_TOP_SRCDIR) \
LALSUITE_DATA_PATH=.:$${LALPULSAR_PREFIX:+"$${LALPULSAR_PREFIX}/share/lalpulsar"}:$${LALSUITE_TOP_SRCDIR:+"$${LALSUITE_TOP_SRCDIR}/lalpulsar/test/"}

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
