EXTRA_DIST = \
	LALVCSInfo.h.in \
	generate_vcs_info.py

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
liblal_la_LIBADD = $(liblal_la_DEPENDENCIES) @LAL_LIBADD@
liblal_la_LDFLAGS = -version-info $(LIBVERSION)

BUILT_SOURCES = \
	._vcsinfo \
	LALVCSInfo.h

bin_PROGRAMS = lal-version
lal_version_SOURCES = lal-version.c
lal_version_LDADD = liblal.la

CLEANFILES =
MOSTLYCLEANFILES = LALBuildDate.c

TESTS = lal-version

LALBuildDate.c: $(liblal_la_DEPENDENCIES)
	$(AM_V_GEN)echo 'const char *lalBuildDate = "'`date +"%Y-%m-%dT%H:%M:%S%z"`'";' >$@

LALVCSInfo.h: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py
	@$(PYTHON) $(srcdir)/generate_vcs_info.py LAL $(srcdir) $$PWD

all-local: include-link
include-link: $(pkginclude_HEADERS)
@INCLUDEMKFILE@

swig_headers = $(pkginclude_HEADERS)
include $(top_srcdir)/swig/swiglal_header.am
