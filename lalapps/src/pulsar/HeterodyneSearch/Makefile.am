BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalapps_pulsar_test.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalapps_ssbtodetector \
	lalapps_SplInter \
	$(END_OF_LINE)

if LALFRAME
bin_PROGRAMS += \
	lalapps_heterodyne_pulsar \
	lalapps_create_pulsar_signal_frame \
	$(END_OF_LINE)
endif

if LALINFERENCE
bin_PROGRAMS += \
	lalapps_pulsar_parameter_estimation_nested \
	lalapps_pulsar_frequency_evolution \
	$(END_OF_LINE)
endif

lalapps_ssbtodetector_SOURCES = ssbtodetector.c

lalapps_pulsar_frequency_evolution_SOURCES = \
	pulsar_frequency_evolution.c \
	ppe_models.c \
	ppe_models.h \
	ppe_utils.c \
	ppe_utils.h \
	$(END_OF_LINE)

lalapps_heterodyne_pulsar_SOURCES = \
	heterodyne_pulsar.c \
	heterodyne_pulsar.h \
	$(END_OF_LINE)

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
	pulsar_parameter_estimation_nested.h \
	$(END_OF_LINE)

lalapps_create_pulsar_signal_frame_SOURCES = \
	create_pulsar_signal_frame.c \
	$(END_OF_LINE)

lalapps_SplInter_SOURCES = \
	SpectralInterpolation.c \
	SpectralInterpolation.h \
	$(END_OF_LINE)

if HAVE_PYTHON
pybin_scripts = \
	lalapps_knope \
	lalapps_knope_result_page \
	lalapps_knope_collate_results \
	lalapps_knope_automation_script \
	$(END_OF_LINE)
pkgpython_PYTHON = \
	pulsarhtmlutils.py \
	pulsarpputils.py \
	knope_utils.py \
	$(END_OF_LINE)
endif

# Add shell test scripts to this variable
test_scripts += test_heterodyne_pulsar.sh

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files += make_frame_cache.py

# These tests require LALFrame
if !LALFRAME
skip_tests += test_heterodyne_pulsar.sh
endif
