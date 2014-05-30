BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

pkgpython_PYTHON = \
	__init__.py \
	git_version.py

EXTRA_DIST = \
	git_version.py

vcs_info_sources = git_version.py

#pkgpython_LTLIBRARIES = \
#	cs_gamma.la
#cs_gamma_la_SOURCES = cs_gamma.c
#cs_gamma_la_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) $(PYTHON_CPPFLAGS)
#cs_gamma_la_CFLAGS = $(AM_CFLAGS) $(CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS)
#cs_gamma_la_LIBADD = $(top_builddir)/src/liblalburst.la
#cs_gamma_la_LDFLAGS = $(AM_LDFLAGS) $(LIBS) $(GSL_LIBS) $(LAL_LIBS) $(PYTHON_LDFLAGS) -module -avoid-version
