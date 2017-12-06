BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	EPSearch.h \
	GenerateBurst.h \
	LALBurstConfig.h \
	LALBurstVCSInfo.h \
	LALBurstVCSInfoHeader.h \
	LIGOLwXMLBurstRead.h \
	SnglBurstUtils.h \
	cs_cosmo.h \
	cs_lambda_cosmo.h

lib_LTLIBRARIES = liblalburst.la
liblalburst_la_SOURCES = \
	EPFilters.c \
	EPSearch.c \
	GenerateBurst.c \
	LALBurstVCSInfo.c \
	LIGOLwXMLBurstRead.c \
	SnglBurstUtils.c \
	cs_lambda_cosmo.c
nodist_liblalburst_la_SOURCES = LALBurstBuildInfo.c

liblalburst_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

DISTCLEANFILES = config.h

vcs_info_sources = LALBurstVCSInfoHeader.h
build_info_source = LALBurstBuildInfo.c
$(build_info_source): $(liblalburst_la_DEPENDENCIES)

EXTRA_DIST = \
	LALBurstConfig.h.in \
	TFR.dox \
	config.h.in \
	$(END_OF_LIST)
