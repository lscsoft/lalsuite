BUILT_SOURCES =
# BUILT_SOURCES = bayestar_cosmology.h
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am

# bayestar_cosmology.h: bayestar_cosmology.py
# 	$(AM_V_GEN)$(PYTHON) $< > $@

if HAVE_CHEALPIX
BAYESTARHDR = bayestar_sky_map.h bayestar_distance.h bayestar_moc.h
BAYESTARSRC = bayestar_sky_map.c bayestar_distance.c bayestar_moc.c
endif

ANALYTICFILES = CorrelationMatrix.csv unimodalMeans.csv bimodalMeans.csv

pkgdata_DATA = $(ANALYTICFILES)

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
	cubic_interp.h \
	$(BAYESTARHDR) \
	$(XMLHDR)

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
	logaddexp.h \
	cubic_interp.c \
	$(BAYESTARSRC) \
	$(XMLSRC)

nodist_liblalinference_la_SOURCES = \
	LALInferenceBuildInfoHeader.h \
	LALInferenceVCSInfo.c \
	$(END_OF_LIST)

noinst_HEADERS = six.h bayestar_cosmology.h omp_interruptible.h

liblalinference_la_CFLAGS = $(AM_CFLAGS) $(HDF5_CFLAGS)
liblalinference_la_CPPFLAGS = $(AM_CPPFLAGS) $(HDF5_CPPFLAGS)
liblalinference_la_LDFLAGS = $(AM_LDFLAGS) $(HDF5_LDFLAGS) $(HDF5_LIBS) -version-info $(LIBVERSION)

LDADD = liblalinference.la

bin_PROGRAMS = \
	lalinference_nest \
	lalinference_burst \
	lalinference_datadump \
	lalinference_bench \
	lalinference_version \
	$(END_OF_LIST)

lalinference_nest_SOURCES = LALInferenceNest.c
lalinference_burst_SOURCES = LALInferenceBurst.c
lalinference_datadump_SOURCES = LALInferenceDataDump.c
lalinference_bench_SOURCES = LALInferenceBench.c
lalinference_version_SOURCES = version.c

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	bayestar_cosmology.py \
	$(ANALYTICFILES) \
	config.h.in

DISTCLEANFILES = config.h

vcs_build_info_source = LALInferenceVCSInfo.c
vcs_info_sources = LALInferenceVCSInfoHeader.h
build_info_source = LALInferenceBuildInfoHeader.h
build_info_dependencies = $(liblalinference_la_DEPENDENCIES)

TESTS = \
	lalinference_version \
	$(END_OF_LIST)
