# Makefile for TDS_isolated

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALFRAME
frame_progs = \
	lalapps_heterodyne_pulsar \
	lalapps_pulsar_parameter_estimation \
	lalapps_create_pulsar_signal_frame

TESTS = test_heterodyne_pulsar.sh

if LALINFERENCE
if LALXML
inference_progs = \
	lalapps_pulsar_parameter_estimation_nested

endif
endif
endif
bin_PROGRAMS = $(frame_progs) $(inference_progs) lalapps_CrabTNHeterodyne

lalapps_heterodyne_pulsar_SOURCES = \
	heterodyne_pulsar.c \
	heterodyne_pulsar.h

lalapps_pulsar_parameter_estimation_SOURCES = \
	pulsar_parameter_estimation.c \
	pulsar_parameter_estimation.h

lalapps_pulsar_parameter_estimation_nested_SOURCES = \
	pulsar_parameter_estimation_nested.c \
	pulsar_parameter_estimation_nested.h \
	ppe_models.c \
	ppe_models.h \
	ppe_likelihood.c \
	ppe_likelihood.h \
	ppe_testing.c \
	ppe_testing.h

lalapps_create_pulsar_signal_frame_SOURCES = \
	create_pulsar_signal_frame.c

lalapps_CrabTNHeterodyne_SOURCES = \
	CrabTNHeterodyne.c \
	HeterodyneCrabPulsar.c \
	HeterodyneCrabPulsar.h

if HAVE_PYTHON
pkgpython_DATA = \
	pulsarpputils.py \
	pulsar_pipeline_utils.py
endif

EXTRA_DIST = \
	test_heterodyne_pulsar.sh \
	make_frame_cache H1response.txt \
	H-CW_Injection-875206560-120.gwf \
	pulsarpputils.py \
	pulsar_pipeline_utils.py

TESTS_ENVIRONMENT = LALPULSAR_DATADIR="$(LALPULSAR_DATADIR)"

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common
