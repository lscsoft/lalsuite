BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

SUBDIRS = \
	bayestar \
	imrtgr \
	io \
	plot \
	tiger \
	$(END_OF_LIST)

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	lalinference_pipe_utils.py \
	nest2pos.py \
	wrapper.py \
	bayespputils.py \
	git_version.py \
	$(END_OF_LIST)

if SWIG_BUILD_PYTHON
pymodule_LTLIBRARIES = _bayespputils.la

_bayespputils_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/lib
_bayespputils_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
_bayespputils_la_LDFLAGS = $(top_builddir)/lib/liblalinference.la -shared -module -avoid-version

endif

pymodule-all-local: _bayespputils.so

_bayespputils.so:
	rm -f $@ && $(LN_S) .libs/$@

CLEANFILES = _bayespputils.so

endif
