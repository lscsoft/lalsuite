BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalapps.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

bin_PROGRAMS = lalapps_ssbtodetector

if LALFRAME
bin_PROGRAMS += \
	lalapps_heterodyne_pulsar \
	lalapps_create_pulsar_signal_frame \
	lalapps_SplInter \
	$(END_OF_LINE)
endif

if LALFRAME
if LALINFERENCE
bin_PROGRAMS += \
	lalapps_pulsar_parameter_estimation_nested \
	lalapps_pulsar_frequency_evolution \
	$(END_OF_LINE)
endif
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
	$(END_OF_LINE)
pkgpython_PYTHON = \
	pulsarhtmlutils.py \
	pulsarpputils.py \
	knope_utils.py \
	$(END_OF_LINE)
endif

EXTRA_DIST += \
	test_heterodyne_pulsar.sh \
	make_frame_cache H1response.txt \
	H-CW_Injection-875206560-120.gwf
