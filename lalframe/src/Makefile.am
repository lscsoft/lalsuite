AM_CPPFLAGS = -I$(top_builddir)/include
LDADD = liblalframe.la

bin_PROGRAMS = \
	lalfr-lld \
	lalfr-cksum \
	lalfr-stat \
	lalfr-dump \
	lalfr-print \
	lalfr-fmt \
	lalfr-cat \
	lalfr-split \
	lalfr-cut \
	lalfr-paste

lalfr_lld_SOURCES = lld.c
lalfr_cksum_SOURCES = cksum.c
lalfr_stat_SOURCES = stat.c
lalfr_dump_SOURCES = dump.c
lalfr_print_SOURCES = print.c
lalfr_fmt_SOURCES = fmt.c
lalfr_cat_SOURCES = cat.c utils.c utils.h
lalfr_split_SOURCES = split.c utils.c utils.h
lalfr_cut_SOURCES = cut.c utils.c utils.h
lalfr_paste_SOURCES = paste.c utils.c utils.h

lalframeincludedir = $(includedir)/lal
lalframeinclude_HEADERS = \
	Aggregation.h \
	FrameCache.h \
	FrameCalibration.h \
	FrameData.h \
	FrameStream.h \
	LALFrameConfig.h \
	LALFrameL.h \
	LALFrameIO.h \
	LALFrameLibVCSInfo.h \
	LALFrameU.h \
	LALFrameVCSInfo.h \
	LowLatencyData.h \
	lalframeswig_test.i

lib_LTLIBRARIES = liblalframe.la
liblalframe_la_SOURCES = \
	LALFrameU.c \
	Aggregation.c \
	FrameCache.c \
	FrameCalibration.c \
	FrameData.c \
	FrameSeries.c \
	FrameStream.c \
	LALFrameIO.c \
	LALFrameVCSInfo.c \
	LowLatencyData.c

noinst_HEADERS = \
	FrameSeriesRead_source.c \
	FrameSeriesWrite_source.c

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
	@echo "# generic user environment script" >$@
	@echo "set LALFRAME_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@printf "$(SWIG_USER_ENV)" >>$@

lalframe-user-env.csh lalframe-user-env.sh: lalframe-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	generate_user_env.awk \
	LALFrameC.c \
	LALFrameConfig.h.in \
	LALFrameL.c

if HAVE_GIT_REPO
BUILT_SOURCES += LALFrameVCSInfo.h
MOSTLYCLEANFILES += LALFrameVCSInfo.h
LALFrameVCSInfo.h .PHONY : gen-LALFrameVCSInfo.h
gen-LALFrameVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALFrameVCSInfo.h $(srcdir)/LALFrameVCSInfo.h.git
endif

@LAL_HEADER_LINK@

swig_headers = $(lalframeinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
