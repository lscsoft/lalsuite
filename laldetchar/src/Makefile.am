BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

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
	$(END_OF_LIST)

nodist_liblaldetchar_la_SOURCES = \
	LALDetCharBuildInfoHeader.h \
	LALDetCharVCSInfo.c \
	$(END_OF_LIST)

liblaldetchar_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

dist_bin_SCRIPTS = \
	laldetchar-idq-append-file-to-cache

LDADD = liblaldetchar.la

bin_PROGRAMS = \
	laldetchar_version \
	$(END_OF_LIST)
#	laldetchar-hveto

#laldetchar_hveto_SOURCES = hveto.c
laldetchar_version_SOURCES = version.c

TESTS = \
	laldetchar_version \
	$(END_OF_LIST)

EXTRA_DIST = \
	LALDetCharConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_build_info_source = LALDetCharVCSInfo.c
vcs_info_sources = LALDetCharVCSInfoHeader.h
build_info_source = LALDetCharBuildInfoHeader.h
build_info_dependencies = $(liblaldetchar_la_DEPENDENCIES)
