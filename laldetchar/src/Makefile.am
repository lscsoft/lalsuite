BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

LDADD = liblaldetchar.la

pkginclude_HEADERS = \
	LALDetChar.h \
	LALDetCharConfig.h \
	LALDetCharGlib.h \
	LALDetCharGlibTypes.h \
	LALDetCharVCSInfo.h \
	LALDetCharRange.h \
	LALDetCharHvetoUtils.h \
	LALDetCharHveto.h \
	LALDetCharVCSInfoHeader.h

lib_LTLIBRARIES = liblaldetchar.la
liblaldetchar_la_SOURCES = \
	LALDetCharGlibTypes.c \
	LALDetCharInspiralRange.c \
	LALDetCharHvetoUtils.c \
	LALDetCharHveto.c \
	LALDetCharVCSInfo.c
nodist_liblaldetchar_la_SOURCES = LALDetCharBuildInfo.c

liblaldetchar_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

#bin_PROGRAMS = \
	#laldetchar-hveto

#laldetchar_hveto_SOURCES = hveto.c

EXTRA_DIST = \
	LALDetCharConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_info_sources = LALDetCharVCSInfoHeader.h
build_info_source = LALDetCharBuildInfo.c
$(build_info_source): $(liblaldetchar_la_DEPENDENCIES)
