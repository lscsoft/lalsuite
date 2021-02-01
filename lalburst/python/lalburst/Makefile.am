BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am
include $(top_srcdir)/gnuscripts/lalsuite_vcs_info.am

vcs_info_sources = git_version.py

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)

pymodule_PYTHON = \
	__init__.py \
	binjfind.py \
	bucluster.py \
	burca.py \
	burca_tailor.py \
	calc_likelihood.py \
	cafe.py \
	date.py \
	git_version.py \
	offsetvector.py \
	packing.py \
	power.py \
	SimBurstUtils.py \
	SnglBurstUtils.py \
	snglcluster.py \
	snglcoinc.py \
	stringutils.py \
	timeslides.py \
	$(END_OF_LIST)

if SWIG_BUILD_PYTHON
pkgpython_LTLIBRARIES = \
	cs_gamma.la
cs_gamma_la_SOURCES = cs_gamma.c
cs_gamma_la_CPPFLAGS = $(AM_CPPFLAGS) $(SWIG_PYTHON_CPPFLAGS)
cs_gamma_la_CFLAGS = $(AM_CFLAGS) $(SWIG_PYTHON_CFLAGS) $(GSL_CFLAGS) -Wno-cast-qual
cs_gamma_la_LDFLAGS = $(top_builddir)/lib/liblalburst.la $(AM_LDFLAGS) $(GSL_LIBS) -module -avoid-version
endif

endif
