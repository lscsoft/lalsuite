AM_CPPFLAGS = @AM_CPPFLAGS@ -I$(top_builddir)/include

bayestardir = $(pkgpythondir)/bayestar

bayestar_PYTHON = \
	__init__.py \
	command.py \
	decorator.py \
	filter.py \
	fits.py \
	ligolw.py \
	ligolw_sky_map.py \
	postprocess.py \
	plot.py \
	timing.py

dist_bayestar_DATA = \
	ne_simplified_coastline.json

if HAVE_CHEALPIX
if SWIG_BUILD_PYTHON
bayestar_LTLIBRARIES = sky_map.la
sky_map_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
sky_map_la_CFLAGS = $(PYTHON_CFLAGS) -w
sky_map_la_LDFLAGS = -shared -module -avoid-version
sky_map_la_LIBADD = $(top_builddir)/src/liblalinference.la
endif
endif
