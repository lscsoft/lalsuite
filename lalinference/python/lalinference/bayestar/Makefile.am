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
	distance.py \
	filter.py \
	ligolw.py \
	ligolw_sky_map.py \
	postprocess.py \
	timing.py \
	$(END_OF_LIST)

if HAVE_CHEALPIX
if SWIG_BUILD_PYTHON
pymodule_LTLIBRARIES = sky_map.la _distance.la

sky_map_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
sky_map_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
sky_map_la_LDFLAGS = -shared -module -avoid-version
sky_map_la_LIBADD = $(top_builddir)/src/liblalinference.la

_distance_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS) -I$(top_srcdir)/src
_distance_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) -Wno-error
_distance_la_LDFLAGS = -shared -module -avoid-version
_distance_la_LIBADD = $(top_builddir)/src/liblalinference.la
endif
endif

all-local: sky_map.so _distance.so

sky_map.so:
	rm -f $@ && $(LN_S) .libs/$@

_distance.so:
	rm -f $@ && $(LN_S) .libs/$@

CLEANFILES = sky_map.so _distance.so

endif
