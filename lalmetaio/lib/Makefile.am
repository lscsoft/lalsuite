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
	LIGOLwXMLRead.h \
	LIGOMetadataTables.h \
	LIGOMetadataUtils.h

lib_LTLIBRARIES = liblalmetaio.la

liblalmetaio_la_SOURCES = \
	LIGOLwXML.c \
	LIGOLwXMLArray.c \
	LIGOLwXMLRead.c \
	LIGOMetadataUtils.c \
	process_params.c \
	processtable.c \
	search_summary.c \
	segment.c \
	sim_burst.c \
	sim_inspiral.c \
	sim_ringdown.c \
	sngl_burst.c \
	sngl_inspiral.c \
	sngl_ringdown.c \
	time_slide.c \
	$(END_OF_LIST)

nodist_liblalmetaio_la_SOURCES = \
	LALMetaIOBuildInfoHeader.h \
	LALMetaIOVCSInfo.c \
	$(END_OF_LIST)

liblalmetaio_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	LALMetaIOConfig.h.in \
	config.h.in

DISTCLEANFILES = config.h

vcs_build_info_source = LALMetaIOVCSInfo.c
vcs_info_sources = LALMetaIOVCSInfoHeader.h
build_info_source = LALMetaIOBuildInfoHeader.h
build_info_dependencies = $(liblalmetaio_la_DEPENDENCIES)
