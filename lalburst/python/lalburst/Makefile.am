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
	cafe.py \
	date.py \
	git_version.py \
	ligolw_burca2.py \
	ligolw_burca_tailor.py \
	packing.py \
	SimBurstUtils.py \
	SnglBurstUtils.py \
	snglcluster.py \
	stringutils.py \
	timeslides.py \
	$(END_OF_LIST)

pkgpython_LTLIBRARIES = \
	cs_gamma.la
cs_gamma_la_SOURCES = cs_gamma.c
cs_gamma_la_CPPFLAGS = $(AM_CPPFLAGS) $(CPPFLAGS) $(PYTHON_CPPFLAGS) -I$(top_srcdir)/src
cs_gamma_la_CFLAGS = $(AM_CFLAGS) $(CFLAGS) $(GSL_CFLAGS) $(LAL_CFLAGS) -Wno-cast-qual
cs_gamma_la_LIBADD = $(top_builddir)/src/liblalburst.la
cs_gamma_la_LDFLAGS = $(AM_LDFLAGS) $(LIBS) $(GSL_LIBS) $(LAL_LIBS) $(PYTHON_LIBS) -module -avoid-version

endif
