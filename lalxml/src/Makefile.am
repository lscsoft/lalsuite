BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	LALXML.h \
	LALXMLConfig.h \
	LALXMLVCSInfo.h \
	LALXMLVCSInfoHeader.h \
	LALXMLVOTableCommon.h \
	LALXMLVOTableSerializers.h

lib_LTLIBRARIES = liblalxml.la
liblalxml_la_SOURCES = \
	LALXML.c \
	LALXMLVCSInfo.c \
	LALXMLVOTableCommon.c \
	LALXMLVOTableSerializers.c
nodist_liblalxml_la_SOURCES = LALXMLBuildInfoHeader.h
liblalxml_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

DISTCLEANFILES = \
	config.h

EXTRA_DIST = \
	config.h.in \
	LALXMLConfig.h.in

vcs_build_info_source = LALXMLVCSInfo.c
vcs_info_sources = LALXMLVCSInfoHeader.h
build_info_source = LALXMLBuildInfoHeader.h
build_info_dependencies = $(liblalxml_la_DEPENDENCIES)

LIBADD = liblalxml.la

bin_PROGRAMS = \
	lalxml_version \
	$(END_OF_LIST)

lalxml_version_SOURCES = version.c

TESTS = \
	lalxml_version \
	$(END_OF_LIST)
