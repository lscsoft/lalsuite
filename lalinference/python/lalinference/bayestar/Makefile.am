BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/bayestar

pymodule_PYTHON = \
	__init__.py \
	command.py \
	decorator.py \
	deprecation.py \
	filter.py \
	ligolw.py \
	sky_map.py \
	postprocess.py \
	timing.py \
	$(END_OF_LIST)

if HAVE_CHEALPIX
if SWIG_BUILD_PYTHON
pymodule_LTLIBRARIES = _sky_map.la

_sky_map_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_sky_map_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
_sky_map_la_LDFLAGS = $(top_builddir)/src/liblalinference.la -shared -module -avoid-version
endif
endif

all-local: _sky_map.so

_sky_map.so:
	rm -f $@ && $(LN_S) .libs/$@

CLEANFILES = _sky_map.so

endif
