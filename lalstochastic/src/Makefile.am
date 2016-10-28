BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	AstroOmega.h \
	CoarseGrainFrequencySeries.h \
	LALStochasticConfig.h \
	LALStochasticVCSInfo.h \
	LALStochasticVCSInfoHeader.h \
	LIGOLwXMLStochasticRead.h \
	SimulatePopcorn.h \
	SimulateSB.h \
	StochasticCrossCorrelation.h \
	StochasticMC.h

lib_LTLIBRARIES = liblalstochastic.la

liblalstochastic_la_SOURCES = \
	AstroOmega.c \
	CoarseGrainFrequencySeries.c \
	LIGOLwXMLStochasticRead.c \
	OverlapReductionFunction.c \
	SimulatePopcorn.c \
	SimulateSB.c \
	StochasticCrossCorrelation.c \
	StochasticInverseNoise.c \
	StochasticMC.c \
	StochasticOmegaGW.c \
	StochasticOptimalFilter.c \
	StochasticOptimalFilterNormalization.c \
	ZeroPadAndFFT.c \
	$(END_OF_LIST)

nodist_liblalstochastic_la_SOURCES = \
	LALStochasticBuildInfoHeader.h \
	LALStochasticVCSInfo.c \
	$(END_OF_LIST)

liblalstochastic_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

DISTCLEANFILES = config.h

vcs_build_info_source = LALStochasticVCSInfo.c
vcs_info_sources = LALStochasticVCSInfoHeader.h
build_info_source = LALStochasticBuildInfoHeader.h
build_info_dependencies = $(liblalstochastic_la_DEPENDENCIES)

bin_PROGRAMS = \
	lalstochastic_version \
	$(END_OF_LIST)

lalstochastic_version_SOURCES = version.c
lalstochastic_version_LDADD = liblalstochastic.la

TESTS = \
	lalstochastic_version \
	$(END_OF_LIST)

EXTRA_DIST = \
	LALStochasticConfig.h.in \
	config.h.in \
	stochastic_CrossCorrFlowchart.png \
	utilitiesCoarseGrain.png \
	$(END_OF_LIST)
