EXTRA_DIST = \
	generate_user_env.awk

lib_LTLIBRARIES = liblal.la

pkginclude_HEADERS = \
	LALLibVCSInfo.h \
	LALVCSInfo.h

liblal_la_SOURCES = \
	LALEmpty.c \
	LALVCSInfo.c
nodist_liblal_la_SOURCES = LALBuildDate.c
liblal_la_DEPENDENCIES = \
	$(top_builddir)/packages/std/src/libstd.la \
	$(top_builddir)/packages/tools/src/libtools.la \
	$(top_builddir)/packages/factories/src/libfactories.la \
	$(top_builddir)/packages/vectorops/src/libvectorops.la \
	$(top_builddir)/packages/utilities/src/libutilities.la \
	$(top_builddir)/packages/stats/src/libstats.la \
	$(top_builddir)/packages/inject/src/libinject.la \
	$(top_builddir)/packages/date/src/libdate.la \
	$(top_builddir)/packages/tdfilter/src/libtdfilter.la \
	$(top_builddir)/packages/window/src/libwindow.la \
	$(top_builddir)/packages/fft/src/libfft.la \
	$(top_builddir)/packages/noisemodels/src/libnoisemodels.la
liblal_la_LIBADD = $(liblal_la_DEPENDENCIES)
liblal_la_LDFLAGS = $(AM_LDFLAGS) -version-info $(LIBVERSION)

bin_PROGRAMS = lal-version
lal_version_SOURCES = lal-version.c
lal_version_LDADD = liblal.la

CLEANFILES = \
	lal-user-env.shell \
	lal-user-env.sh \
	lal-user-env.csh

sysconf_DATA = \
	lal-user-env.sh \
	lal-user-env.csh

BUILT_SOURCES = \
	lal-user-env.sh \
	lal-user-env.csh

MOSTLYCLEANFILES = LALBuildDate.c

TESTS = lal-version

LALBuildDate.c: $(liblal_la_DEPENDENCIES)
	$(AM_V_GEN)echo 'const char *lalBuildDate = "'`date +"%Y-%m-%dT%H:%M:%S%z"`'";' >$@

lal-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# source this file to access $(PACKAGE_NAME)" >$@
	@echo "set LAL_PREFIX $(prefix)" >>$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PKG_CONFIG_PATH $(libdir)/pkgconfig" >>$@
	@echo "prepend PYTHONPATH $(pyexecdir) $(pythondir)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@
	@echo "prepend OCTAVE_PATH $(octexecdir)" >>$@

lal-user-env.csh lal-user-env.sh: lal-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED="$(SED)" -v output="$@" $<

if HAVE_GIT_REPO
BUILT_SOURCES += LALVCSInfo.h
MOSTLYCLEANFILES += LALVCSInfo.h
LALVCSInfo.h .PHONY : gen-LALVCSInfo.h
gen-LALVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALVCSInfo.h $(srcdir)/LALVCSInfo.h.git
endif

all-local: include-link
include-link: $(pkginclude_HEADERS)
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
