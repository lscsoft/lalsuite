BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
dist_bin_SCRIPTS =
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_help2man.am
include $(top_srcdir)/gnuscripts/lalsuite_old_exec_name_stubs.am

if FRAMEL
if LALFRAME
FRAMESRCS = \
	series.c \
	series.h \
	lalappsfrutils.c \
	lalappsfrutils.h
endif
endif

noinst_LTLIBRARIES = liblalapps.la

liblalapps_la_SOURCES = \
	LALAppsVCSInfo.h \
	LALAppsVCSInfoHeader.h \
	$(FRAMESRCS)

nodist_liblalapps_la_SOURCES = \
	LALAppsBuildInfoHeader.h \
	LALAppsVCSInfo.c \
	$(END_OF_LIST)

bin_PROGRAMS = lalapps_version

lalapps_version_SOURCES = version2.c

LDADD = liblalapps.la

if HAVE_PYTHON
pkgpython_PYTHON = \
	__init__.py \
	git_version.py
endif

TESTS = lalapps_version

vcs_build_info_source = LALAppsVCSInfo.c
vcs_info_sources = LALAppsVCSInfoHeader.h git_version.py
build_info_source = LALAppsBuildInfoHeader.h
build_info_dependencies = $(liblalapps_la_DEPENDENCIES)

EXTRA_DIST += \
	__init__.py \
	git_version.py \
	lalappsfrutils.c \
	lalappsfrutils.h \
	series.c \
	series.h \
	$(END_OF_LIST)

old_exec_name_stubs = \
	lalapps_cache \
	lalapps_fftw_wisdom \
	lalapps_fftwf_wisdom \
	lalapps_path2cache \
	lalapps_searchsum2cache \
	lalapps_tconvert \
	$(END_OF_LIST)

new_exec_name_regex = 's|^lalapps_|lal_|'
