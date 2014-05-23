BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

if LALXML
XMLHDR = LALInferenceXML.h
XMLSRC = LALInferenceXML.c
endif

if HAVE_CHEALPIX
BAYESTARHDR = bayestar_sky_map.h
BAYESTARSRC = bayestar_sky_map.c
endif

pkginclude_HEADERS = \
	LALInference.h \
	LALInferenceConfig.h \
	LALInferenceLikelihood.h \
	LALInferenceNestedSampler.h \
	LALInferencePrior.h \
	LALInferenceReadData.h \
	LALInferenceTemplate.h \
	LALInferenceProposal.h \
	LALInferenceVCSInfo.h \
	LALInferenceInit.h \
	LALInferenceVCSInfoHeader.h \
	LALInferenceRemoveLines.h \
	$(BAYESTARHDR) \
	$(XMLHDR)

lib_LTLIBRARIES = liblalinference.la
liblalinference_la_SOURCES = \
	LALInference.c \
	LALInferenceLikelihood.c \
	LALInferenceNestedSampler.c \
	LALInferencePrior.c \
	LALInferenceReadData.c \
	LALInferenceTemplate.c \
	LALInferenceProposal.c \
	LALInferenceInitCBC.c \
	LALInferenceVCSInfo.c \
	LALInferenceRemoveLines.c \
	logaddexp.h \
	likely.h \
	$(BAYESTARSRC) \
	$(XMLSRC)
nodist_liblalinference_la_SOURCES = LALInferenceBuildInfo.c

liblalinference_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALInferenceConfig.h.in \
	LALInferenceXML.c \
	LALInferenceXML.h \
	chapter_inference.dox \
	config.h.in

DISTCLEANFILES = config.h

vcs_info_sources = LALInferenceVCSInfoHeader.h
build_info_source = LALInferenceBuildInfo.c
$(build_info_source): $(liblalinference_la_DEPENDENCIES)
