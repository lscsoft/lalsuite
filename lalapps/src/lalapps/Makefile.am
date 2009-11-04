if FRAME
FRAMESRCS = \
	series.c \
	series.h \
	lalappsfrutils.c \
	lalappsfrutils.h
else
FRAMESRCS =
endif

noinst_LTLIBRARIES = liblalapps.la
liblalapps_la_SOURCES = \
	lalapps.c \
	lalapps.h \
	getopt.c \
	getopt1.c \
	getopt.h \
	processtable.c \
	processtable.h \
	$(FRAMESRCS)

bin_PROGRAMS = \
	lalapps_version \
	lalapps_tconvert \
	lalapps_cache

lalapps_version_SOURCES = \
	version.c \
	lalapps.h \
	LALAppsVCSInfo.h
lalapps_tconvert_SOURCES = \
	tconvert.c \
	getdate.c \
	getdate.h
lalapps_cache_SOURCES = \
	cache.c

LDADD = liblalapps.la @LAL_LIBS@
INCLUDES = -I$(top_srcdir)/src

if HAVE_PYTHON
pkgpython_DATA = \
	__init__.py \
	defaults.py
endif

EXTRA_DIST = \
	series.c \
	series.h \
	getdate.y \
	__init__.py \
	defaults.py \
	lalappsfrutils.c \
	lalappsfrutils.h
	git_version \
	LALAppsVCSInfo.h.in

BUILT_SOURCE = LALAppsVCSInfo.h

MOSTLYCLEANFILES = .vcsinfo

.vcsinfo:
	@touch .vcsinfo

LALAppsVCSInfo.h: LALAppsVCSInfo.h.in .vcsinfo git_version
	@echo "creating LALAppsVCSInfo.h"
	@buildd="`pwd`" && cd $(srcdir) && sourced="`pwd`" && cd $$buildd && \
		./git_version --sed LALApps $$sourced $$buildd
	@-rm -f .vcsinfo


TESTS = lalapps_version
