if FRAME
if LALFRAME
FRAMESRCS = \
	series.c \
	series.h \
	lalappsfrutils.c \
	lalappsfrutils.h
endif
endif

if LALMETAIO
LALMETAIOSRCS = \
	processtable.c \
	processtable.h
endif

if FFTW
lalapps_fftwf_wisdom_SOURCES = \
	fftwf_wisdom.c
lalapps_fftw_wisdom_SOURCES = \
	fftw_wisdom.c
FFTWPROGS = \
	lalapps_fftwf_wisdom \
	lalapps_fftw_wisdom
endif

noinst_LTLIBRARIES = liblalapps.la
liblalapps_la_SOURCES = \
	lalapps.c \
	lalapps.h \
	getopt.c \
	getopt1.c \
	getopt.h \
	LALAppsVCSInfo.c \
	LALAppsVCSInfo.h \
	$(FRAMESRCS) \
	$(LALMETAIOSRCS)

bin_PROGRAMS = \
	lalapps_version \
	lalapps_tconvert \
	lalapps_cache \
	$(FFTWPROGS)

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
pkgpython_PYTHON = \
	__init__.py \
	git_version.py
dist_bin_SCRIPTS = \
	lalapps_path2cache
endif

EXTRA_DIST = \
	series.c \
	series.h \
	getdate.y \
	__init__.py \
	fftw_wisdom.c \
	fftwf_wisdom.c \
	lalappsfrutils.c \
	lalappsfrutils.h \
	processtable.c \
	processtable.h \
	LALAppsVCSInfo.h.in \
	generate_vcs_info.py\
	git_version.py.in

BUILT_SOURCES = \
	._vcsinfo \
	LALAppsVCSInfo.h \
	git_version.py

LALAppsVCSInfo.h: ._vcsinfo

git_version.py: ._vcsinfo

._vcsinfo: $(srcdir)/generate_vcs_info.py
	@$(PYTHON) $(srcdir)/generate_vcs_info.py LALApps $(srcdir) $$PWD

TESTS = lalapps_version
