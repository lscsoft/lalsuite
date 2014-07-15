BUILT_SOURCES =
MOSTLYCLEANFILES =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

pymoduledir = $(pkgpythondir)/bayestar

pymodule_PYTHON = \
	__init__.py \
	command.py \
	decorator.py \
	filter.py \
	fits.py \
	ligolw.py \
	ligolw_sky_map.py \
	plot.py \
	postprocess.py \
	timing.py \
	$(END_OF_LIST)

if HAVE_CHEALPIX
if SWIG_BUILD_PYTHON
pymodule_LTLIBRARIES = sky_map.la
sky_map_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS)
sky_map_la_CFLAGS = $(SWIG_PYTHON_CFLAGS) -w
sky_map_la_LDFLAGS = -shared -module -avoid-version
sky_map_la_LIBADD = $(top_builddir)/src/liblalinference.la
endif
endif
