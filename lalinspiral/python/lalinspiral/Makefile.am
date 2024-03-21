BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

.PHONY: thinca-all-local

all-local: thinca-all-local

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	git_version.py \
	inspinjfind.py \
	thinca.py \
	$(END_OF_LIST)

if SWIG_BUILD_PYTHON

pymodule_LTLIBRARIES = _thinca.la

_thinca_la_SOURCES = _thinca.c
_thinca_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -DMODULE_NAME="\"lalinspiral._thinca\""
_thinca_la_CFLAGS = $(AM_CFLAGS) -fno-strict-aliasing -DMODULE_NAME="\"lalinspiral._thinca\""
_thinca_la_LDFLAGS = $(AM_LDFLAGS) -module -avoid-version $(SWIG_PYTHON_FLAT_NAMESPACE_LDFLAG)

thinca-all-local: _thinca.so

_thinca.so:
	rm -f $@ && $(LN_S) .libs/$@

CLEANFILES = _thinca.so

endif

endif
