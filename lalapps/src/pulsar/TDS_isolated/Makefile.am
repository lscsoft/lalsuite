# Makefile for TDS_isolated 
# $Id$

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALFRAME
  bin_PROGRAMS = lalapps_pulsar_parameter_estimation lalapps_heterodyne_pulsar
else
  bin_PROGRAMS =
endif

lalapps_pulsar_parameter_estimation_SOURCES = pulsar_parameter_estimation.c pulsar_parameter_estimation.h

lalapps_heterodyne_pulsar_SOURCES = heterodyne_pulsar.c heterodyne_pulsar.h

EXTRA_DIST = test_heterodyne_pulsar.sh make_frame_cache H1response.txt H-CW_Injection-875206560-120.gwf

TESTS = test_heterodyne_pulsar.sh

TESTS_ENVIRONMENT = \
LAL_TOP_SRCDIR=$(LAL_TOP_SRCDIR) \
LAL_PREFIX=$(LAL_PREFIX) \
LAL_DATA_PATH=.:$${LAL_TOP_SRCDIR:+"$${LAL_TOP_SRCDIR}/packages/pulsar/test"}:$${LAL_PREFIX:+"$${LAL_PREFIX}/share/lal"}
