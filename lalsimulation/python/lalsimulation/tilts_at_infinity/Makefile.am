BUILT_SOURCES =
MOSTLYCLEANFILES =
EXTRA_DIST =
include $(top_srcdir)/gnuscripts/lalsuite_python.am

if HAVE_PYTHON

pymoduledir = $(pkgpythondir)/tilts_at_infinity

pymodule_PYTHON = \
	__init__.py \
	calc_tilts_prec_avg_regularized.py \
	hybrid_spin_evolution.py \
	tilts_at_infinity_utils.py \
	$(END_OF_LIST)

endif
