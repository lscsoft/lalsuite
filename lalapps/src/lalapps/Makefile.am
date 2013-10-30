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
AM_CPPFLAGS = -I$(top_srcdir)/src -DPKG_DATA_DIR='"$(pkgdatadir)/"'

if HAVE_PYTHON
pkgpython_PYTHON = \
	__init__.py \
	git_version.py
dist_bin_SCRIPTS = \
	lalapps_path2cache \
	lalapps_searchsum2cache
endif

EXTRA_DIST = \
	series.c \
	series.h \
	getdate.y \
	__init__.py \
	fftw_wisdom.c \
	fftwf_wisdom.c \
	generate_user_env.awk \
	lalappsfrutils.c \
	lalappsfrutils.h \
	processtable.c \
	processtable.h \
	git_version.py

MOSTLYCLEANFILES =

CLEANFILES = \
	lalapps-user-env.shell \
	lalapps-user-env.sh \
	lalapps-user-env.csh

sysconf_DATA = \
	lalapps-user-env.sh \
	lalapps-user-env.csh

BUILT_SOURCES = \
	lalapps-user-env.sh \
	lalapps-user-env.csh

if HAVE_PYTHON
lalapps_pythonpath = $(pyexecdir) $(pythondir)
else
lalapps_pythonpath =
endif

lalapps-user-env.shell: Makefile $(srcdir)/generate_user_env.awk
	@echo "# source this file to access $(PACKAGE_NAME)" >$@
	@echo "prepend PATH $(bindir)" >>$@
	@echo "prepend MANPATH $(mandir)" >>$@
	@echo "prepend PYTHONPATH $(lalapps_pythonpath)" >>$@
	@echo "prepend LAL_DATA_PATH $(pkgdatadir)" >>$@

lalapps-user-env.csh lalapps-user-env.sh: lalapps-user-env.shell $(srcdir)/generate_user_env.awk
	$(AM_V_GEN)$(AWK) -f $(srcdir)/generate_user_env.awk -v SED="$(SED)" -v output="$@" $<

if HAVE_GIT_REPO
BUILT_SOURCES += LALAppsVCSInfo.h
MOSTLYCLEANFILES += LALAppsVCSInfo.h
LALAppsVCSInfo.h .PHONY : gen-LALAppsVCSInfo.h
gen-LALAppsVCSInfo.h :
	$(AM_V_at)$(GENERATE_VCS_INFO) LALAppsVCSInfo.h $(srcdir)/LALAppsVCSInfo.h.git
endif

if HAVE_GIT_REPO
BUILT_SOURCES += git_version.py
MOSTLYCLEANFILES += git_version.py
git_version.py .PHONY : gen-git_version.py
gen-git_version.py :
	$(AM_V_at)$(GENERATE_VCS_INFO) git_version.py $(srcdir)/git_version.py.git
endif

TESTS = lalapps_version
