BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lallibs_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	CLR.h \
	EPSearch.h \
	GenerateBurst.h \
	LALBurstConfig.h \
	LALBurstVCSInfo.h \
	LALBurstVCSInfoHeader.h \
	LIGOLwXMLBurstRead.h \
	LIGOMetadataBurstUtils.h \
	TFTransform.h \
	TSData.h \
	TSSearch.h \
	Thresholds.h \
	TimeFreq.h \
	TrackSearch.h \
	cs_cosmo.h \
	cs_lambda_cosmo.h \
	lalburstswig_test.i

lib_LTLIBRARIES = liblalburst.la
liblalburst_la_SOURCES = \
	CleanAll.c \
	CreateTFPlane.c \
	CreateTimeFreqParam.c \
	CreateTimeFreqRep.c \
	DestroyTimeFreqParam.c \
	DestroyTimeFreqRep.c \
	Dwindow.c \
	EPSearch.c \
	FreqSeriesToTFPlane.c \
	GenerateBurst.c \
	HarmonicFinder.c \
	LALBurstVCSInfo.c \
	LIGOLwXMLBurstRead.c \
	RefInterference.c \
	SnglBurstUtils.c \
	TSData.c \
	TfrPswv.c \
	TfrRsp.c \
	TfrSp.c \
	TfrWv.c \
	Thresholds.c \
	TrackSearch.c \
	cs_lambda_cosmo.c
nodist_liblalburst_la_SOURCES = LALBurstBuildInfo.c

liblalburst_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALBurstConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_info_sources = LALBurstVCSInfoHeader.h
build_info_source = LALBurstBuildInfo.c
$(build_info_source): $(liblalburst_la_DEPENDENCIES)
