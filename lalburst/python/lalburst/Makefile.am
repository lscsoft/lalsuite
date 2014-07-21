BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

vcs_info_sources = git_version.py

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	git_version.py \
	$(END_OF_LIST)

#pkgpython_LTLIBRARIES = \
#	cs_gamma.la
#cs_gamma_la_SOURCES = cs_gamma.c
#cs_gamma_la_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) $(PYTHON_CPPFLAGS)
#cs_gamma_la_CFLAGS = $(AM_CFLAGS) $(CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS)
#cs_gamma_la_LIBADD = $(top_builddir)/src/liblalburst.la
#cs_gamma_la_LDFLAGS = $(AM_LDFLAGS) $(LIBS) $(GSL_LIBS) $(LAL_LIBS) $(PYTHON_LDFLAGS) -module -avoid-version
