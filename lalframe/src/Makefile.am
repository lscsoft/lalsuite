INCLUDES = -I$(top_builddir)/include
LDADD = liblalframe.la

bin_PROGRAMS = \
	lalframe_lld

lalframe_lld_SOURCES = lld.c

lalframeincludedir = $(includedir)/lal
lalframeinclude_HEADERS = \
	Aggregation.h \
	FrameCache.h \
	FrameCalibration.h \
	FrameData.h \
	FrameStream.h \
	LALFrameConfig.h \
	LALFrameIO.h \
	LALFrameL.h \
	LALFrameLibVCSInfo.h \
	LALFrameVCSInfo.h \
	LowLatencyData.h

lib_LTLIBRARIES = liblalframe.la
liblalframe_la_SOURCES = \
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

CLEANFILES = \
	lalframe-user-env.shell \
	lalframe-user-env.sh \
	lalframe-user-env.csh

sysconf_DATA = \
	lalframe-user-env.sh \
	lalframe-user-env.csh

BUILT_SOURCES = \
	._vcsinfo \
	LALFrameVCSInfo.h \
	lalframe-user-env.sh \
	lalframe-user-env.csh

lalframe-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# generic user environment script" >$@
	@echo "set LALFRAME_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend LD_LIBRARY_PATH $(libdir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@printf "$(SWIG_USER_ENV)" >>$@

lalframe-user-env.csh lalframe-user-env.sh: lalframe-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED=$(SED) -v package=$(PACKAGE_NAME) -v output=$@ $< || { rm -f $@; exit 1; }

EXTRA_DIST = \
	generate_user_env.awk \
	generate_vcs_info.py \
	LALFrameConfig.h.in \
	LALFrameVCSInfo.h.in

LALFrameVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py LALFrameConfig.h
	@$(PYTHON) $(srcdir)/generate_vcs_info.py LALFrame $(srcdir) $$PWD

@LAL_HEADER_LINK@

swig_headers = $(lalframeinclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
