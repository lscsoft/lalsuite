BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lallibs_header_links.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

lib_LTLIBRARIES = liblal.la

pkginclude_HEADERS = \
	LALConfig.h \
	LALVCSInfo.h \
	LALVCSInfoHeader.h

liblal_la_SOURCES = \
	LALEmpty.c \
	LALVCSInfo.c
nodist_liblal_la_SOURCES = LALBuildInfo.c
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

TESTS = lal-version

vcs_info_sources = LALVCSInfoHeader.h
build_info_source = LALBuildInfo.c
$(build_info_source): $(liblal_la_DEPENDENCIES)
