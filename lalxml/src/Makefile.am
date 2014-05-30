BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lallibs_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	LALXML.h \
	LALXMLConfig.h \
	LALXMLVCSInfo.h \
	LALXMLVCSInfoHeader.h \
	LALXMLVOTableCommon.h \
	LALXMLVOTableSerializers.h \
	lalxmlswig_test.i

lib_LTLIBRARIES = liblalxml.la
liblalxml_la_SOURCES = \
	LALXML.c \
	LALXMLVCSInfo.c \
	LALXMLVOTableCommon.c \
	LALXMLVOTableSerializers.c
nodist_liblalxml_la_SOURCES = LALXMLBuildInfo.c
liblalxml_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

DISTCLEANFILES = \
	config.h

EXTRA_DIST = \
	config.h.in \
	LALXMLConfig.h.in

vcs_info_sources = LALXMLVCSInfoHeader.h
build_info_source = LALXMLBuildInfo.c
$(build_info_source): $(liblalxml_la_DEPENDENCIES)
