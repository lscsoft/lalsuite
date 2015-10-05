BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/imrtgr

pymodule_PYTHON = \
	__init__.py \
	imrtgrutils.py \
	nrutils.py \
	pneqns.py \
	$(END_OF_LIST)

webpagedir = $(pkgdatadir)/imrtgr_webpage_templates

webpage_DATA = \
	result.html \
	webpage.css \
	webpage.js \
	jquery.js \
	$(END_OF_LIST)

priordatadir = $(pkgdatadir)/imrtgr_prior_data

priordata_DATA = \
	Prior_Mfchif_nonprecspin_Healy2014_comp_mass_min1.0_comp_mass_max500.0_comp_spin_min-1.0_comp_spin_max1.0.pklz \
	$(END_OF_LIST)
	nrutils.py \
	$(END_OF_LIST)

endif
