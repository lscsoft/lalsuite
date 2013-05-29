AM_CPPFLAGS = -I$(top_builddir)/include

pkgpython_PYTHON = \
	__init__.py \
	git_version.py

#pkgpython_LTLIBRARIES = \
#	cs_gamma.la

EXTRA_DIST = \
	git_version.py

BUILT_SOURCES =
MOSTLYCLEANFILES =

if HAVE_GIT_REPO
BUILT_SOURCES += git_version.py
MOSTLYCLEANFILES += git_version.py
git_version.py .PHONY : gen-git_version.py
gen-git_version.py :
	$(AM_V_at)$(GENERATE_VCS_INFO) git_version.py $(srcdir)/git_version.py.git
endif

#cs_gamma_la_SOURCES = cs_gamma.c
#cs_gamma_la_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) $(PYTHON_CPPFLAGS)
#cs_gamma_la_CFLAGS = $(AM_CFLAGS) $(CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS)
#cs_gamma_la_LIBADD = $(top_builddir)/src/liblalburst.la
#cs_gamma_la_LDFLAGS = $(AM_LDFLAGS) $(LIBS) $(GSL_LIBS) $(LAL_LIBS) $(PYTHON_LDFLAGS) -module -avoid-version
