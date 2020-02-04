BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	sbank \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	git_version.py \
	inspinjfind.py \
	thinca.py \
	$(END_OF_LIST)

pymodule_LTLIBRARIES = _thinca.la

_thinca_la_SOURCES = _thinca.c
_thinca_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -DMODULE_NAME="\"lalinspiral._thinca\""
_thinca_la_CFLAGS = $(AM_CFLAGS) -fno-strict-aliasing -DMODULE_NAME="\"lalinspiral._thinca\""
_thinca_la_LDFLAGS = $(AM_LDFLAGS) $(SWIG_PYTHON_LIBS) -module -avoid-version

pymodule-all-local: _thinca.so

_thinca.so:
	rm -f $@ && $(LN_S) .libs/$@

CLEANFILES = _thinca.so

endif
