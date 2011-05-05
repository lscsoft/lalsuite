# Makefile for TDS_isolated

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALFRAME
bin_PROGRAMS = \
	lalapps_heterodyne_pulsar \
	lalapps_pulsar_parameter_estimation 

TESTS = test_heterodyne_pulsar.sh
endif

lalapps_heterodyne_pulsar_SOURCES = \
	heterodyne_pulsar.c \
	heterodyne_pulsar.h

lalapps_pulsar_parameter_estimation_SOURCES = \
	pulsar_parameter_estimation.c \
	pulsar_parameter_estimation.h

#lalapps_pulsar_parameter_estimation_nested_SOURCES = \
#	parameter_estimation_nested.c \
#	parameter_estimation_nested.h \
#	LALInference.c  LALInference.h \
#	LALInferenceNestedSampler.c  LALInferenceNestedSampler.h

EXTRA_DIST = \
	test_heterodyne_pulsar.sh \
	make_frame_cache H1response.txt \
	H-CW_Injection-875206560-120.gwf

TESTS_ENVIRONMENT = \
LALSUITE_TOP_SRCDIR=$(LALSUITE_TOP_SRCDIR) \
LALSUITE_DATA_PATH=.:$${LALPULSAR_PREFIX:+"$${LALPULSAR_PREFIX}/share/lalpulsar"}:$${LALSUITE_TOP_SRCDIR:+"$${LALSUITE_TOP_SRCDIR}/lalpulsar/test/"}

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
