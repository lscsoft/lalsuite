BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	AstroOmega.h \
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
	LALStochasticVCSInfo.c \
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
	ZeroPadAndFFT.c
nodist_liblalstochastic_la_SOURCES = LALStochasticBuildInfo.c
liblalstochastic_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALStochasticConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_info_sources = LALStochasticVCSInfoHeader.h
build_info_source = LALStochasticBuildInfo.c
$(build_info_source): $(liblalstochastic_la_DEPENDENCIES)
