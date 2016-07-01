include $(top_srcdir)/gnuscripts/lalapps.am
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

# Makefile for HeterodyneSearch

if LALFRAME
frame_progs = \
	lalapps_heterodyne_pulsar \
	lalapps_pulsar_parameter_estimation \
	lalapps_create_pulsar_signal_frame \
	lalapps_SplInter

TESTS = test_heterodyne_pulsar.sh

if LALINFERENCE
inference_progs = \
	lalapps_pulsar_parameter_estimation_nested \
	lalapps_pulsar_frequency_evolution
endif
endif
bin_PROGRAMS = $(frame_progs) $(inference_progs) lalapps_CrabTNHeterodyne lalapps_ssbtodetector

lalapps_ssbtodetector_SOURCES = ssbtodetector.c

lalapps_pulsar_frequency_evolution_SOURCES = \
	pulsar_frequency_evolution.c \
	ppe_models.c \
	ppe_models.h \
	ppe_utils.c \
	ppe_utils.h

lalapps_heterodyne_pulsar_SOURCES = \
	heterodyne_pulsar.c \
	heterodyne_pulsar.h

lalapps_pulsar_parameter_estimation_SOURCES = \
	pulsar_parameter_estimation.c \
	pulsar_parameter_estimation.h

lalapps_pulsar_parameter_estimation_nested_SOURCES = \
	ppe_models.c \
	ppe_models.h \
	ppe_likelihood.c \
	ppe_likelihood.h \
	ppe_testing.c \
	ppe_testing.h \
	ppe_init.c \
	ppe_init.h \
	ppe_utils.c \
	ppe_utils.h \
	ppe_readdata.c \
	ppe_readdata.h \
	ppe_inject.c \
	ppe_inject.h \
	ppe_roq.c \
	ppe_roq.h \
	pulsar_parameter_estimation_nested.c \
	pulsar_parameter_estimation_nested.h

lalapps_create_pulsar_signal_frame_SOURCES = \
	create_pulsar_signal_frame.c

lalapps_CrabTNHeterodyne_SOURCES = \
	CrabTNHeterodyne.c \
	HeterodyneCrabPulsar.c \
	HeterodyneCrabPulsar.h

lalapps_SplInter_SOURCES = \
	SpectralInterpolation.c \
	SpectralInterpolation.h

if HAVE_PYTHON
pybin_scripts = \
	lalapps_knope \
	lalapps_knope_create_prior \
	lalapps_knope_result_page \
	lalapps_knope_collate_results
pkgpython_PYTHON = \
	pulsarhtmlutils.py \
	pulsarpputils.py \
	knope_utils.py
endif

EXTRA_DIST += \
	test_heterodyne_pulsar.sh \
	make_frame_cache H1response.txt \
	H-CW_Injection-875206560-120.gwf
