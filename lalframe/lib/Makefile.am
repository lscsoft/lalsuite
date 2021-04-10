BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

FRAMECSRCS = LALFrameUFrameC.h LALFrameC.c
FRAMELSRCS = LALFrameUFrameL.h LALFrameL.c
FRAMEUSRCS =

if FRAMEC
FRAMEUSRCS += $(FRAMECSRCS)
endif
if FRAMEL
FRAMEUSRCS += $(FRAMELSRCS)
endif

LDADD = liblalframe.la

pkginclude_HEADERS = \
	FrameCalibration.h \
	LALFrStream.h \
	LALFrameConfig.h \
	LALFrameIO.h \
	LALFrameU.h \
	LALFrameVCSInfo.h \
	LALFrameVCSInfoHeader.h

lib_LTLIBRARIES = liblalframe.la

liblalframe_la_SOURCES = \
	$(FRAMEUSRCS) \
	LALFrameU.c \
	LALFrameIO.c \
	LALFrStream.c \
	LALFrStreamRead.c \
	LALFrStreamLegacy.c \
	FrameCalibration.c \
	$(END_OF_LIST)

nodist_liblalframe_la_SOURCES = \
	LALFrameBuildInfoHeader.h \
	LALFrameVCSInfo.c \
	$(END_OF_LIST)

noinst_HEADERS = \
	LALFrStreamReadTS_source.c \
	LALFrStreamReadFS_source.c \
	LALFrameIO_source.c

liblalframe_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

EXTRA_DIST = \
	$(FRAMECSRCS) \
	$(FRAMELSRCS) \
	config.h.in \
	LALFrameC.c \
	LALFrameL.c \
	LALFrameConfig.h.in

DISTCLEANFILES = config.h

vcs_build_info_source = LALFrameVCSInfo.c
vcs_info_sources = LALFrameVCSInfoHeader.h
build_info_source = LALFrameBuildInfoHeader.h
build_info_dependencies = $(liblalframe_la_DEPENDENCIES)
