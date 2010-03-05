EXTRA_DIST = \
	LALVCSInfo.h.in \
	git_version

lib_LTLIBRARIES = liblal.la

pkginclude_HEADERS = LALVCSInfo.h

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

BUILT_SOURCES = LALVCSInfo.h

bin_PROGRAMS = lal-version
lal_version_SOURCES = lal-version.c
lal_version_LDADD = liblal.la

MOSTLYCLEANFILES = \
	.vcsinfo \
	LALBuildDate.c \
	git_info.sed

TESTS = lal-version

LALBuildDate.c: $(liblal_la_DEPENDENCIES)
	@echo 'const char *lalBuildDate = "'`date +"%Y-%m-%dT%H:%M:%S%z"`'";' >$@

.vcsinfo:
	@touch .vcsinfo

LALVCSInfo.h: LALVCSInfo.h.in .vcsinfo git_version
	@echo "creating LALVCSInfo.h"
	@buildd="`pwd`" && cd $(srcdir) && sourced="`pwd`" && cd $$buildd && \
	./git_version --sed LAL $$sourced $$buildd
	@-rm -f .vcsinfo

all-local: include-link
include-link: LALVCSInfo.h
	@d=`pwd`; \
		targetdir=`cd $(top_builddir)/include/lal && pwd`; \
		if test ! -r $$targetdir/LALVCSInfo.h ; then \
			rm -f $$targetdir/LALVCSInfo.h; \
			$(LN_S) $$d/LALVCSInfo.h $$targetdir; \
		fi
