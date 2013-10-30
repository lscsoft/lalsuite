AM_CPPFLAGS = -I$(top_builddir)/include -DPKG_DATA_DIR='"$(pkgdatadir)/"'
LDADD = liblalframe.la

bin_PROGRAMS = \
	lalfr-cksum \
	lalfr-stat \
	lalfr-dump \
	lalfr-print \
	lalfr-fmt \
	lalfr-cat \
	lalfr-split \
	lalfr-cut \
	lalfr-paste \
	lalfr-stream \
	lalfr-vis
man_MANS = \
	lalfr-cksum.1 \
	lalfr-stat.1 \
	lalfr-dump.1 \
	lalfr-print.1 \
	lalfr-fmt.1 \
	lalfr-cat.1 \
	lalfr-split.1 \
	lalfr-cut.1 \
	lalfr-paste.1 \
	lalfr-stream.1 \
	lalfr-vis.1

lalfr_cksum_SOURCES = cksum.c
lalfr_stat_SOURCES = stat.c
lalfr_dump_SOURCES = dump.c
lalfr_print_SOURCES = print.c
lalfr_fmt_SOURCES = fmt.c
lalfr_cat_SOURCES = cat.c utils.c utils.h
lalfr_split_SOURCES = split.c utils.c utils.h
lalfr_cut_SOURCES = cut.c utils.c utils.h
lalfr_paste_SOURCES = paste.c utils.c utils.h
lalfr_stream_SOURCES = stream.c
lalfr_vis_SOURCES = vis.c

lalframeincludedir = $(includedir)/lal
lalframeinclude_HEADERS = \
	Aggregation.h \
	FrameCalibration.h \
	LALFrStream.h \
	LALFrameConfig.h \
	LALFrameIO.h \
	LALFrameL.h \
	LALFrameU.h \
	LALFrameLibVCSInfo.h \
	LALFrameVCSInfo.h \
	lalframeswig_test.i

lib_LTLIBRARIES = liblalframe.la
liblalframe_la_SOURCES = \
	LALFrameU.c \
	LALFrameIO.c \
	LALFrStream.c \
	LALFrStreamRead.c \
	LALFrStreamLegacy.c \
	LALFrameVCSInfo.c \
	FrameCalibration.c \
	Aggregation.c

noinst_HEADERS = \
	LALFrStreamReadTS_source.c \
	LALFrStreamReadFS_source.c \
	LALFrameIO_source.c

liblalframe_la_LDFLAGS = -version-info $(LIBVERSION)

MOSTLYCLEANFILES = \
	lalframe-user-env.shell \
	lalframe-user-env.sh \
	lalframe-user-env.csh

sysconf_DATA = \
	lalframe-user-env.sh \
	lalframe-user-env.csh

BUILT_SOURCES = \
	lalframe-user-env.sh \
	lalframe-user-env.csh

lalframe-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# source this file to access $(PACKAGE_NAME)" >$@
	@echo "set LALFRAME_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lalframe-user-env.csh lalframe-user-env.sh: lalframe-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED="$(SED)" -v output="$@" $<

EXTRA_DIST = \
	$(MANS) \
	config.h.in \
	generate_user_env.awk \
	LALFrameC.c \
	LALFrameL.c \
	LALFrameConfig.h.in

if HAVE_GIT_REPO
BUILT_SOURCES += LALFrameVCSInfo.h
MOSTLYCLEANFILES += LALFrameVCSInfo.h
LALFrameVCSInfo.h .PHONY : gen-LALFrameVCSInfo.h
gen-LALFrameVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALFrameVCSInfo.h $(srcdir)/LALFrameVCSInfo.h.git
endif

DISTCLEANFILES = config.h

@LAL_HEADER_LINK@

swig_headers = $(lalframeinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
