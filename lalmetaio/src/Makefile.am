BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkginclude_HEADERS = \
	LALMetaIOConfig.h \
	LALMetaIOVCSInfo.h \
	LALMetaIOVCSInfoHeader.h \
	LIGOLwXML.h \
	LIGOLwXMLArray.h \
	LIGOLwXMLHeaders.h \
	LIGOLwXMLInspiralHeaders.h \
	LIGOLwXMLRead.h \
	LIGOMetadataTables.h \
	LIGOMetadataUtils.h

lib_LTLIBRARIES = liblalmetaio.la
liblalmetaio_la_SOURCES = \
	CreateMetaTableDir.c \
	LALMetaIOVCSInfo.c \
	LIGOLwXML.c \
	LIGOLwXMLArray.c \
	LIGOLwXMLRead.c \
	LIGOMetadataUtils.c
nodist_liblalmetaio_la_SOURCES = LALMetaIOBuildInfoHeader.h

liblalmetaio_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALMetaIOConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_build_info_source = LALMetaIOVCSInfo.c
vcs_info_sources = LALMetaIOVCSInfoHeader.h
build_info_source = LALMetaIOBuildInfoHeader.h
build_info_dependencies = $(liblalmetaio_la_DEPENDENCIES)

LIBADD = liblalmetaio.la

bin_PROGRAMS = \
	lalmetaio_version \
	$(END_OF_LIST)

lalmetaio_version_SOURCES = version.c

TESTS = \
	lalmetaio_version \
	$(END_OF_LIST)
