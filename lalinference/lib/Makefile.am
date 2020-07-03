BUILT_SOURCES =
# BUILT_SOURCES = bayestar_cosmology.h
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

# bayestar_cosmology.h: bayestar_cosmology.py
#	$(AM_V_GEN)$(PYTHON) $< > $@

ANALYTICFILES = \
	CorrelationMatrix.csv \
	bimodalMeans.csv \
	unimodalMeans.csv

FIDUCIALINJECTIONS = \
	fiducialBBH.xml \
	fiducialBNS.xml

EXAMPLE_CONFIG = \
	config_tiger_example.ini \
	lalinference_pipe_example.ini \
	lib_pipe_example.ini

pkgdata_DATA = \
	$(ANALYTICFILES) \
	$(EXAMPLE_CONFIG) \
	$(FIDUCIALINJECTIONS)

pkginclude_HEADERS = \
	LALInference.h \
	LALInferenceConfig.h \
	LALInferenceLikelihood.h \
	LALInferenceMultibanding.h \
	LALInferenceNestedSampler.h \
	LALInferencePrior.h \
	LALInferenceReadBurstData.h \
	LALInferenceReadData.h \
	LALInferenceTemplate.h \
	LALInferenceProposal.h \
	LALInferenceClusteredKDE.h \
	LALInferenceKDE.h \
	LALInferenceVCSInfo.h \
	LALInferenceInit.h \
	LALInferenceVCSInfoHeader.h \
	LALInferenceRemoveLines.h \
	LALInferenceCalibrationErrors.h \
	LALInferenceGenerateROQ.h \
	LALInferenceBurstRoutines.h \
	LALInferenceHDF5.h \
	LALInferencePriorVolumes.h \
	LALInferenceDistanceMarg.h \
	cubic_interp.h \
	distance_integrator.h

lib_LTLIBRARIES = liblalinference.la

liblalinference_la_SOURCES = \
	LALInference.c \
	LALInferenceLikelihood.c \
	LALInferenceAnalyticLikelihood.c \
	LALInferenceMultibanding.c \
	LALInferenceNestedSampler.c \
	LALInferencePrior.c \
	LALInferenceReadBurstData.c \
	LALInferenceReadData.c \
	LALInferenceTemplate.c \
	LALInferenceProposal.c \
	LALInferenceClusteredKDE.c \
	LALInferenceKDE.c \
	LALInferenceInitBurst.c \
	LALInferenceInitCBC.c \
	LALInferenceRemoveLines.c \
	LALInferenceCalibrationErrors.c \
	LALInferenceGenerateROQ.c \
	LALInferenceBurstRoutines.c \
	LALInferenceHDF5.c \
	LALInferencePriorVolumes.c \
	DetectorFixedSkyCoords.c \
	LALInferenceDistanceMarg.c \
	logaddexp.h \
	cubic_interp.c \
	distance_integrator.c

nodist_liblalinference_la_SOURCES = \
	LALInferenceBuildInfoHeader.h \
	LALInferenceVCSInfo.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	bayestar_cosmology.h \
	omp_interruptible.h \
	six.h

liblalinference_la_CFLAGS = $(AM_CFLAGS) $(HDF5_CFLAGS)
liblalinference_la_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)
liblalinference_la_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LDFLAGS) $(HDF5_LIBS) -version-info $(LIBVERSION)

LDADD = liblalinference.la

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	bayestar_cosmology.py \
	$(ANALYTICFILES) \
	$(FIDUCIALINJECTIONS) \
	$(EXAMPLE_CONFIG) \
	config.h.in

DISTCLEANFILES = config.h

vcs_build_info_source = LALInferenceVCSInfo.c
vcs_info_sources = LALInferenceVCSInfoHeader.h
build_info_source = LALInferenceBuildInfoHeader.h
build_info_dependencies = $(liblalinference_la_DEPENDENCIES)
