BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

LDADD = liblalinference.la

lalinference_nest_SOURCES = LALInferenceNest.c
lalinference_burst_SOURCES = LALInferenceBurst.c
lalinference_datadump_SOURCES = LALInferenceDataDump.c
lalinference_bench_SOURCES = LALInferenceBench.c

if HAVE_CHEALPIX
BAYESTARHDR = bayestar_sky_map.h bayestar_distance.h
BAYESTARSRC = bayestar_sky_map.c bayestar_distance.c
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
	LALInferenceVCSInfo.c \
	LALInferenceRemoveLines.c \
	LALInferenceCalibrationErrors.c \
	LALInferenceGenerateROQ.c \
	LALInferenceBurstRoutines.c \
	LALInferenceHDF5.c \
	LALInferencePriorVolumes.c \
	DetectorFixedSkyCoords.c \
	logaddexp.h \
	$(BAYESTARSRC) \
	$(XMLSRC)
nodist_liblalinference_la_SOURCES = LALInferenceBuildInfo.c

liblalinference_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

bin_PROGRAMS = \
	lalinference_nest \
	lalinference_burst \
	lalinference_datadump \
	lalinference_bench

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	$(ANALYTICFILES) \
	config.h.in

DISTCLEANFILES = config.h

vcs_info_sources = LALInferenceVCSInfoHeader.h
build_info_source = LALInferenceBuildInfo.c
$(build_info_source): $(liblalinference_la_DEPENDENCIES)
