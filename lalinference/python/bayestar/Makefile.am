AM_CPPFLAGS = -I$(top_builddir)/include

# The Numpy 1.7 headers contain some static function definitions that
# themselves are not '-Werror'-safe. Turn off the '-Werror' (treat warnings
# as errors) flag, for just targets declared in this file.
CFLAGS += -Wno-error

if HAVE_PYTHON
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

bayestar_DATA = \
	ne_simplified_coastline.json

if HAVE_CHEALPIX
if SWIG_BUILD_PYTHON
bayestar_LTLIBRARIES = sky_map.la
sky_map_la_CPPFLAGS = $(AM_CPPFLAGS) $(PYTHON_CPPFLAGS)
sky_map_la_CFLAGS = $(AM_CFLAGS) $(PYTHON_CFLAGS)
sky_map_la_LDFLAGS = $(AM_LDFLAGS) -shared -module -avoid-version
sky_map_la_LIBADD = $(LALSUPPORT_LIBS) $(top_builddir)/src/liblalinference.la
endif
endif

endif
