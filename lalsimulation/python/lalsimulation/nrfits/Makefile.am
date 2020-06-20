BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/nrfits

pymodule_PYTHON = \
	__init__.py \
	eval_fits.py \
	nrfits.py \
	NRSur7dq4Remnant.py \
	NRSur3dq8Remnant.py \
	pn_spin_evolution_wrapper.py \
	quaternion_utils.py \
	$(END_OF_LIST)

endif
