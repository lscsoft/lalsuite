BUILT_SOURCES =
EXTRA_DIST =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalpulsar_bin.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am

bin_PROGRAMS = \
	lalpulsar_ssbtodetector \
	$(END_OF_LINE)

noinst_PROGRAMS = \
	lalpulsar_SplInter \
	$(END_OF_LINE)

if LALFRAME

bin_PROGRAMS += \
	lalpulsar_heterodyne \
	$(END_OF_LINE)

noinst_PROGRAMS += \
	lalpulsar_create_signal_frame \
	$(END_OF_LINE)

endif

if LALSIMULATION
if LALINFERENCE

bin_PROGRAMS += \
	lalpulsar_frequency_evolution \
	lalpulsar_parameter_estimation_nested \
	$(END_OF_LINE)

endif
endif

lalpulsar_ssbtodetector_SOURCES = \
	ssbtodetector.c \
	$(END_OF_LINE)

lalpulsar_frequency_evolution_SOURCES = \
	ppe_models.c \
	ppe_models.h \
	ppe_utils.c \
	ppe_utils.h \
	pulsar_frequency_evolution.c \
	$(END_OF_LINE)

lalpulsar_heterodyne_SOURCES = \
	heterodyne_pulsar.c \
	heterodyne_pulsar.h \
	$(END_OF_LINE)

lalpulsar_parameter_estimation_nested_SOURCES = \
	ppe_init.c \
	ppe_init.h \
	ppe_inject.c \
	ppe_inject.h \
	ppe_likelihood.c \
	ppe_likelihood.h \
	ppe_models.c \
	ppe_models.h \
	ppe_readdata.c \
	ppe_readdata.h \
	ppe_roq.c \
	ppe_roq.h \
	ppe_testing.c \
	ppe_testing.h \
	ppe_utils.c \
	ppe_utils.h \
	pulsar_parameter_estimation_nested.c \
	pulsar_parameter_estimation_nested.h \
	$(END_OF_LINE)

lalpulsar_create_signal_frame_SOURCES = \
	create_pulsar_signal_frame.c \
	$(END_OF_LINE)

lalpulsar_SplInter_SOURCES = \
	SpectralInterpolation.c \
	SpectralInterpolation.h \
	$(END_OF_LINE)

if HAVE_PYTHON

pybin_scripts = \
	lalpulsar_knope \
	lalpulsar_knope_automation_script \
	lalpulsar_knope_collate_results \
	lalpulsar_knope_result_page \
	$(END_OF_LINE)

pkgpython_PYTHON = \
	knope_utils.py \
	pulsarhtmlutils.py \
	pulsarpputils.py \
	$(END_OF_LINE)

endif

# Add shell test scripts to this variable
test_scripts += test_frequency_evolution.sh
if HAVE_PYTHON
test_scripts += test_heterodyne.sh
test_scripts += test_ssbtodetector.py
test_scripts += test_flat_prior.py
test_scripts += test_flat_evidence_default.py
#test_scripts += test_flat_evidence_walk.py
endif

# Add any helper programs required by tests to this variable
test_helpers +=

# Add any extra files required by tests (e.g. helper scripts) to this variable
test_extra_files += make_frame_cache.py

# These tests require LALFrame
if !LALFRAME
skip_tests += test_heterodyne.sh
endif

# These tests require LALInference
if !LALINFERENCE
skip_tests += test_frequency_evolution.sh
skip_tests += test_flat_prior.py
skip_tests += test_flat_evidence_default.py
skip_tests += test_flat_evidence_walk.py
endif
