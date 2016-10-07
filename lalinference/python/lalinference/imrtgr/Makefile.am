BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

FINALBHFITS = \
	Prior_Mfaf_nonprec_Healy2014_M1-500_aligned.pklz \
	Prior_Mfaf_nonprec_Healy2014_M1-500_isotropic.pklz \
	$(END_OF_LIST)

WEBPAGE = \
	result.html \
	webpage.css \
	webpage.js \
	jquery.js \
	$(END_OF_LIST)

EXTRA_DIST += \
	$(FINALBHFITS) \
	$(WEBPAGE) \
	$(END_OF_LIST)

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/imrtgr

pymodule_PYTHON = \
	__init__.py \
	imrtgrutils.py \
	nrutils.py \
	pneqns.py \
	$(END_OF_LIST)

webpagedir = $(pkgdatadir)/imrtgr_webpage_templates
webpage_DATA = $(WEBPAGE)

priordatadir = $(pkgdatadir)/imrtgr_prior_data
priordata_DATA = $(FINALBHFITS)

endif
